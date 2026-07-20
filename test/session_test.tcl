mol new receptor.pdbqt type pdb waitfor all
mol rename top "receptor"
mol modstyle 0 0 NewCartoon
mol modcolor 0 0 ColorID 8
mol modmaterial 0 0 Transparent

mol new biased/biased_out.pdbqt type pdb waitfor all
mol rename top "biased_poses"
mol modstyle 0 1 Licorice 0.3 12 12
mol modcolor 0 1 ColorID 10

mol new conventional/conventional_out.pdbqt type pdb waitfor all
mol rename top "conventional_poses"
mol modstyle 0 2 Licorice 0.3 12 12
mol modcolor 0 2 ColorID 3

draw color yellow
draw sphere {1.782 66.634 6.337} radius 0.5 resolution 20

color Display Background white
display projection Orthographic
display resetview
puts "OK_ZALADOWANO: [molinfo num] molekul"
