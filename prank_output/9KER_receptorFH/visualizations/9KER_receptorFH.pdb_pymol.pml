from pymol import cmd,stored

set depth_cue, 1
set fog_start, 0.4

set_color b_col, [36,36,85]
set_color t_col, [10,10,10]
set bg_rgb_bottom, b_col
set bg_rgb_top, t_col      
set bg_gradient

set  spec_power  =  200
set  spec_refl   =  0

load "data/9KER_receptorFH.pdb", protein
create ligands, protein and organic
select xlig, protein and organic
delete xlig

hide everything, all

color white, elem c
color bluewhite, protein
#show_as cartoon, protein
show surface, protein
#set transparency, 0.15

show sticks, ligands
set stick_color, magenta




# SAS points

load "data/9KER_receptorFH.pdb_points.pdb.gz", points
hide nonbonded, points
show nb_spheres, points
set sphere_scale, 0.2, points
cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)


stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
lastSTP=stored.list[-1] # get the index of the last residue
hide lines, resn STP

cmd.select("rest", "resn STP and resi 0")

for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))



set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [3270,3271,3272,3276,3277,3280,3648,3649,3654,3655,3650,5178,5192,5193,5194,4780,5204,4797,4799,4800,4809,3934,3935,3936,4487,4489,4491,4760,3929,4774,3695,3697,3698,4779,4785,4393,4395,5691,5695,5696,5697,5698,5154,5155,5167,3265,5164,3354,3353,5990,5680,5685,5146,5991,5992,5996,5147,5241,5203,5198,5199,4781,5220,4798,5242,3341,5997,6028,6031] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.616,0.278,0.702]
select surf_pocket2, protein and id [5678,5684,5989,5990,5685,5687,5721,5991,5996,994,996,1027,1037,1039,1901,416,444,2547,2545,2546,382,383,384,5960,5969,5971,5983,995,2556,623,2583,2584,2599,2600,5997,2557] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.902,0.361,0.361]
select surf_pocket3, protein and id [3591,3628,4300,4144,3429,4122] 
set surface_color,  pcol3, surf_pocket3 
   

deselect

orient
