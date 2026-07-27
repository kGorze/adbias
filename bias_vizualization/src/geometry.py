import math

from .models import AutoDockGrid, Bias, BiasGeometry, Edge3D, Index3D, Point3D, SampledPoint


#typy
type AxisRange = tuple[int, int]
type GridRanges = tuple[AxisRange, AxisRange, AxisRange]

#pomocniczne funkcje
def _clamp_int(value: int, minimum: int, maximum: int) -> int:
    """
    ogranicza wartość całkowitą do przedziału [minimum, maximum]"""
    return max(minimum, min(value, maximum))

def distance3(first: Point3D, second: Point3D) -> float:
    """
    dystans euklidesowaski z math.dist() między dwoma punktami w przestrzeni 3D
    """
    return math.dist(first, second)

def add3(first: Point3D, second: Point3D) -> Point3D:
    """
    dodawanie punktu w przestrzeni 3D, przez ich wartości
    """
    values_of_point = zip(first, second, strict=True)
    added_values = tuple(left + right for left, right in values_of_point)
    return added_values


def nearest_grid_index(
    coordinate: float,
    minimum: float,
    spacing: float,
    point_count: int,
) -> int:
    """
    z punktów centrum, minimum, odstepu siadki i liczby punktów siatki, zwraca indeks najbliższego punktu siatki. ale tylko w jednej osi!
    """

    #trzeba uważać, bo math.floor() zaokrągla w dół
    how_far_from_start = coordinate - minimum
    position_in_grid = how_far_from_start / spacing + 0.5

    unbounded_index = math.floor(position_in_grid)
    index = _clamp_int(unbounded_index, 0, point_count - 1)

    return index

#rysowanie krawędzi pudełka
def _box_edges(minimum: Point3D, maximum: Point3D) -> tuple[Edge3D, ...]:
    """
    zwraca 12 krawędzi w zafixowanym na osi pudełku 3D
    """
    x0, y0, z0 = minimum
    x1, y1, z1 = maximum

    # nie robimy dowolnej ilości tych punktów, dokładniejsza informacja dla type checkera [Point3D, ...]
    first_side: tuple[Point3D, Point3D, Point3D, Point3D] = (
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0)
    )

    second_side: tuple[Point3D, Point3D, Point3D, Point3D] = (
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1)
    )

    vertices: tuple[
        Point3D, Point3D, Point3D, Point3D,
        Point3D, Point3D, Point3D, Point3D
    ] = first_side + second_side

    edge_indices = (
        (0, 1), (1, 2), (2, 3), (3, 0),  # dolna ściana
        (4, 5), (5, 6), (6, 7), (7, 4),  # górna ściana
        (0, 4), (1, 5), (2, 6), (3, 7),  # krawędzie pionowe
    )

    return tuple(
        (vertices[start], vertices[end])
        for start, end in edge_indices
    )


#zmiana na dokładne wyrażenia matematyczne, żeby widzieć co się dzieje
def _corrected_axis_range(
    bias_coordinate: float,
    epsilon_radius: float,
    grid_minimum: float,
    spacing: float,
    point_count: int,
) -> tuple[int, int] | None:
    """
    zwraca indeksy pierwszego i ostatniego punktu siatki wzdłuż osi, który mieści się w sferze epsilon_radius wokół bias_coordinate. Robi to z tolerancą, która jest tutaj lokalną zmienną!
    """

    tolerance = 1.0e-12

    x_min = grid_minimum
    c_minus_re = bias_coordinate - epsilon_radius
    c_plus_re = bias_coordinate + epsilon_radius
    h = spacing

    # Punkt siatki o indeksie i ma współrzędną:
    #     x_i = x_min + i*h
    #
    # Szukamy wszystkich punktów należących do przedziału:
    #     C - R_epsilon <= x_i <= C + R_epsilon
    #
    # Po podstawieniu x_i i rozwiązaniu względem indeksu i:
    #     (C - R_epsilon - x_min) / h
    #         <= i <=
    #     (C + R_epsilon - x_min) / h
    #
    # Ponieważ i musi być liczbą całkowitą:
    #     first = ceil((C - R_epsilon - x_min) / h)
    #     last  = floor((C + R_epsilon - x_min) / h)

    first = math.ceil(
        (c_minus_re - x_min) / h - tolerance
    )
    last = math.floor(
        (c_plus_re - x_min) / h + tolerance
    )

    # Ograniczamy indeksy do rzeczywistego zakresu siatki:
    #     0 <= i <= point_count - 1
    first = max(first, 0)
    last = min(last, point_count - 1)

    if first > last:
        return None

    return first, last


def _nearest_grid_point(
    grid: AutoDockGrid,
    position: Point3D,
) -> tuple[Index3D, Point3D]:
    """Zwraca indeks i współrzędne najbliższego punktu siatki."""

    #używamy już poprzedniej funkcji jaka została zrobiona
    indices = (
        nearest_grid_index(
            position[0],
            grid.minimum[0],
            grid.spacing,
            grid.point_counts[0],
        ),
        nearest_grid_index(
            position[1],
            grid.minimum[1],
            grid.spacing,
            grid.point_counts[1],
        ),
        nearest_grid_index(
            position[2],
            grid.minimum[2],
            grid.spacing,
            grid.point_counts[2],
        ),
    )

    #punkt z trzech wymiarów
    point = (
        grid.minimum[0] + indices[0] * grid.spacing,
        grid.minimum[1] + indices[1] * grid.spacing,
        grid.minimum[2] + indices[2] * grid.spacing,
    )

    return indices, point

def _gaussian_geometry(
    bias: Bias,
    nearest_point: Point3D,
    epsilon: float,
) -> tuple[float, float, float]:
    """
    Zwraca geometrię gaussowskiego biasu:
    odległość środka biasu od najbliższego punktu siatki,
    wartość biasu w tym punkcie oraz promień R_epsilon,
    przy którym bezwzględna wartość biasu spada do epsilon.
    """

    # Oznaczenia matematyczne:
    #     C      - środek biasu
    #     P      - najbliższy punkt siatki
    #     d      - odległość |C - P|
    #     V_set  - wartość biasu w jego centrum
    #     r      - promień / parametr szerokości Gaussa
    #     epsilon - próg, poniżej którego bias uznajemy za pomijalny
    #     R_epsilon - odległość, dla której |V(R_epsilon)| = epsilon

    C = bias.center
    P = nearest_point
    V_set = bias.vset
    r = bias.radius

    # Odległość środka biasu od najbliższego punktu siatki:
    #     d = |C - P|
    # czyli w 3D:
    #     d = sqrt(
    #         (C_x - P_x)^2
    #       + (C_y - P_y)^2
    #       + (C_z - P_z)^2
    #     )
    d = distance3(C, P)

    # Gaussowski bias ma postać:
    #     V(d) = V_set * e^(-d^2 / r^2)
    # Dlatego wartość biasu w najbliższym punkcie siatki wynosi:
    #     V_nearest = V_set * exp(-d^2 / r^2)
    d_squared = d**2
    r_squared = r**2

    V_nearest = V_set * math.exp(
        -d_squared / r_squared
    )
    # Szukamy promienia R_epsilon, dla którego:
    #     |V(R_epsilon)| = epsilon
    # Z równania Gaussa:
    #     epsilon
    #         = |V_set| * exp(-R_epsilon^2 / r^2)
    # otrzymujemy:
    #     epsilon / |V_set|
    #         = exp(-R_epsilon^2 / r^2)
    #
    #     ln(|V_set| / epsilon)
    #         = R_epsilon^2 / r^2
    # a więc:
    #     R_epsilon
    #         = r * sqrt(ln(|V_set| / epsilon))
    R_epsilon = r * math.sqrt(
        math.log(abs(V_set) / epsilon)
    )

    return d, V_nearest, R_epsilon

def _current_grid_box(
    grid: AutoDockGrid,
    nearest_indices: Index3D,
    nearest_point: Point3D,
    radius: float,
) -> tuple[int, float, tuple[Edge3D, ...], GridRanges]:

    #liczba przedziałów od środka do krawędzi
    half_intervals = math.ceil(
        2.0 * radius / grid.spacing
    )

    #rzeczywisty półbok szceścianu
    half_side = half_intervals * grid.spacing

    #minimum i maximum sześcianu
    minimum = tuple(
        value - half_side
        for value in nearest_point
    )
    maximum = tuple(
        value + half_side
        for value in nearest_point
    )


    edges = _box_edges(minimum, maximum)
    ranges = tuple(
        (
            _clamp_int(
                index - half_intervals,
                0,
                point_count - 1,
            ),
            _clamp_int(
                index + half_intervals,
                0,
                point_count - 1,
            ),
        )
        for index, point_count in zip(
            nearest_indices,
            grid.point_counts,
            strict=True,
        )
    )

    return half_intervals, half_side, edges, ranges

def _corrected_grid_box(
    grid: AutoDockGrid,
    center: Point3D,
    epsilon_radius: float,
) -> tuple[Edge3D, ...]:
    """
    Zwraca krawędzie najmniejszego boxa z punktów siatki,
    które mieszczą się w zakresie epsilon_radius wokół center.
    """
    C = center
    R_epsilon = epsilon_radius
    h = grid.spacing

    # Dla każdej osi wyznaczamy indeksy:
    #
    #     i_min <= i <= i_max
    #
    # odpowiadające punktom siatki należącym do zakresu
    #     C - R_epsilon <= x_i <= C + R_epsilon
    ranges = tuple(
        _corrected_axis_range(
            C_axis,
            R_epsilon,
            x_min,
            h,
            point_count,
        )
        for C_axis, x_min, point_count in zip(
            C,
            grid.minimum,
            grid.point_counts,
            strict=True,
        )
    )

    if any(axis_range is None for axis_range in ranges):
        return ()

    complete_ranges = tuple(
        axis_range
        for axis_range in ranges
        if axis_range is not None
    )

    # Punkt siatki o indeksie i ma współrzędną:
    #
    #     x_i = x_min + i*h
    #
    # więc granice boxa wynikają z pierwszego i ostatniego indeksu.
    B_min = tuple(
        x_min + i_min * h
        for x_min, (i_min, _) in zip(
            grid.minimum,
            complete_ranges,
            strict=True,
        )
    )
    B_max = tuple(
        x_min + i_max * h
        for x_min, (_, i_max) in zip(
            grid.minimum,
            complete_ranges,
            strict=True,
        )
    )

    return _box_edges(B_min, B_max)


def _sample_bias_points(
    grid: AutoDockGrid,
    bias: Bias,
    ranges: GridRanges,
    epsilon: float,
) -> tuple[tuple[SampledPoint, ...], int, float]:
    """
    Oblicza wartość gaussowskiego biasu w punktach siatki
    i zaznacza punkty, dla których ΔV < -epsilon.
    """
    x_range, y_range, z_range = ranges

    bx, by, bz = bias.center
    V_set = bias.vset
    r = bias.radius
    h = grid.spacing

    x_min, y_min, z_min = grid.minimum

    points: list[SampledPoint] = []
    accepted_count = 0
    largest_accepted_distance = 0.0

    r_squared = r**2

    for k in range(z_range[0], z_range[1] + 1):
        z = z_min + k * h
        for j in range(y_range[0], y_range[1] + 1):
            y = y_min + j * h
            for i in range(x_range[0], x_range[1] + 1):
                x = x_min + i * h
                # d² = (x-bx)² + (y-by)² + (z-bz)²
                d_squared = (
                    (x - bx) ** 2
                    + (y - by) ** 2
                    + (z - bz) ** 2
                )
                d = math.sqrt(d_squared)
                # f(d) = exp(-d²/r²)
                fraction = math.exp(
                    -d_squared / r_squared
                )
                # ΔV(d) = V_set * f(d)
                delta_V = V_set * fraction
                # Dla atrakcyjnego biasu:
                #     ΔV < -epsilon
                accepted = delta_V < -epsilon

                points.append(
                    SampledPoint(
                        position=(x, y, z),
                        distance=d,
                        fraction=fraction,
                        energy_delta=delta_V,
                        accepted=accepted,
                    )
                )

                if accepted:
                    accepted_count += 1
                    largest_accepted_distance = max(
                        largest_accepted_distance,
                        d,
                    )

    return (
        tuple(points),
        accepted_count,
        largest_accepted_distance,
    ) 

def calculate_bias_geometry(
    grid: AutoDockGrid,
    bias: Bias,
    epsilon: float = 0.01,
) -> BiasGeometry:
    if not math.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("epsilon must be a finite positive number")

    amplitude = abs(bias.vset)
    if amplitude <= epsilon:
        raise ValueError("abs(Vset) must be greater than epsilon")

    point_counts = grid.point_counts
    grid_minimum = grid.minimum

    # poprawka na bardziej czytelny kod w nearest_grid_point()
    nearest_indices, nearest_point = _nearest_grid_point(
    grid,
    bias.center,
    )

    #nowy punkt geometrii gausowskiej
    (center_grid_distance, energy_at_nearest_point, epsilon_radius,
    ) = _gaussian_geometry(bias, nearest_point, epsilon,)

    #liczenie obecnego pudełka siatki
    (half_intervals, half_side, current_box_edges, current_ranges,) = _current_grid_box( grid, nearest_indices, nearest_point, bias.radius,)

    #poprawne pudełko zależne od epsilona
    corrected_box_edges = _corrected_grid_box(
    grid,
    bias.center,
    epsilon_radius,
    )

    #samplowanie punktów
    (
        candidate_points,
        accepted_count,
        largest_accepted_distance,
    ) = _sample_bias_points(
        grid,
        bias,
        current_ranges,
        epsilon,
    )

    return BiasGeometry(
        grid=grid,
        bias=bias,
        epsilon=epsilon,
        nearest_indices=nearest_indices,
        nearest_point=nearest_point,
        center_grid_distance=center_grid_distance,
        energy_at_nearest_point=energy_at_nearest_point,
        epsilon_radius=epsilon_radius,
        current_half_intervals=half_intervals,
        current_half_side=half_side,
        current_box_edges=current_box_edges,
        corrected_box_edges=corrected_box_edges,
        candidate_points=tuple(candidate_points),
        accepted_count=accepted_count,
        largest_accepted_distance=largest_accepted_distance,
    )
