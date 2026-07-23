# === receptor (statyczny) ===
mol new receptor.pdbqt type pdb waitfor all
mol rename top "receptor"
mol modstyle 0 0 NewCartoon
mol modcolor 0 0 ColorID 8
mol modmaterial 0 0 Transparent
# jesli NewCartoon wyswietli sie jako gladka rurka zamiast helis/wstazek,
# to STRIDE (przypisywanie struktury drugorzedowej) nie zadzialal w Twoim
# srodowisku - zamien wtedy powyzsza linie modstyle na: mol modstyle 0 0 Tube

# === 9 niezaleznych poz BIASED jako "trajektoria" (klatki) ===
mol new biased/biased_out.pdbqt type pdb waitfor all
mol rename top "biased_poses"
mol modstyle 0 1 Licorice 0.3 12 12
mol modcolor 0 1 ColorID 10

# === 9 niezaleznych poz CONVENTIONAL jako druga "trajektoria" ===
mol new conventional/conventional_out.pdbqt type pdb waitfor all
mol rename top "conventional_poses"
mol modstyle 0 2 Licorice 0.3 12 12
mol modcolor 0 2 ColorID 3

# === znacznik miejsca biasu (zolta kulka) ===
draw color yellow
draw sphere {1.782 66.634 6.337} radius 0.5 resolution 20

# === wyglad sceny ===
color Display Background white
display projection Orthographic
display resetview
