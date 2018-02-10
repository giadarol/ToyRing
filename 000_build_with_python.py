import numpy as np


L_halfcell = 50.
phase_adv_cell = np.pi/3
n_cells_arc = 16 # to have zero dispersionin the SS needs to be a multiple of 4 
n_arcs = 4
n_dip_half_cell = 3
n_regcells_straight = 6
n_cells_insertion = 6

frac_q_x = .27
frac_q_y = .295
Qpx = 15.
Qpy = 20.

Qx = n_arcs*(n_cells_arc+n_regcells_straight+4+n_cells_insertion)*phase_adv_cell/2./np.pi + frac_q_x
Qy = n_arcs*(n_cells_arc+n_regcells_straight+4+n_cells_insertion)*phase_adv_cell/2./np.pi + frac_q_y

#Assumes DS made by an empty cell and a full cell
circum = n_arcs*(n_cells_arc+n_regcells_straight+4+n_cells_insertion)*L_halfcell*2.
n_dip_total = n_arcs*(n_cells_arc+2)*n_dip_half_cell*2

# Bending angle
b_ang = 2.*np.pi/n_dip_total

# Strength quadrupole
focal_length = L_halfcell/np.sin(phase_adv_cell*1.111111) # I add 10% to avoid integer resonance
k1l_quad = 1./focal_length

# Start building sequence

sequence = ''

sequence+='''
!Quad strengths (to be matched)
kqf:=%e;'''%k1l_quad
sequence+='''
kqd:=%e;'''%(-k1l_quad)

sequence+='''

!Sext strengths (to be matched)
ksf:=0.01;''' # match needs to start from a non zero value
sequence+='''
ksd:=0.01;''' # match needs to start from a non zero value

sequence+='''

mb: multipole,knl=%e;'''%b_ang
sequence+='''
qf: multipole,knl:={0,kqf};
qd: multipole,knl:={0,kqd};
sf: multipole,knl:={0,0,ksf};
sd: multipole,knl:={0,0,ksd};'''


sequence+='''
kqf1:=kqf;
kqd1:=kqd;
kqf2:=kqf;
kqd2:=kqd;
kqf3:=kqf;
kqd3:=kqd;
kqf4:=kqf;
kqd4:=kqd;
qf1: multipole,knl:={0,kqf1};
qd1: multipole,knl:={0,kqd1};
qf2: multipole,knl:={0,kqf2};
qd2: multipole,knl:={0,kqd2};
qf3: multipole,knl:={0,kqf3};
qd3: multipole,knl:={0,kqd3};
qf4: multipole,knl:={0,kqf4};
qd4: multipole,knl:={0,kqd4};'''

sequence+='''

circum = %e;
toyring: sequence, refer=centre, l=circum; 

'''%circum

s_start_cell = 0.
for i_arc in xrange(n_arcs):


	#Dispesion suppressor - full cell
	sequence += 'qf_arc%d_dsleftf: qf, at=%e;\n'%(i_arc, s_start_cell)	
	for i_bend in xrange(n_dip_half_cell):
		sequence += 'mb_arc%d_dsleft_%d: mb, at=%e;\n'%(i_arc, i_bend, 
								s_start_cell+(i_bend+1)*L_halfcell/(n_dip_half_cell+1))
	sequence += 'qd_arc%d_dsleftf: qd, at=%e;\n'%(i_arc, s_start_cell+L_halfcell)
	for i_bend in xrange(n_dip_half_cell):
		sequence += 'mb_arc%d_dsleft_%d: mb, at=%e;\n'%(i_arc, i_bend+n_dip_half_cell, 
								s_start_cell+(i_bend+1)*L_halfcell/(n_dip_half_cell+1)+L_halfcell)
	s_start_cell += L_halfcell*2
	#Dispesion suppressor - empty cell
	sequence += 'qf_arc%d_dslefte: qf, at=%e;\n'%(i_arc, s_start_cell)	
	sequence += 'qd_arc%d_dslefte: qd, at=%e;\n'%(i_arc, s_start_cell+L_halfcell)
	s_start_cell += L_halfcell*2
	
	#Arc
	for i_cell in xrange(n_cells_arc):
		
		sequence += 'qf_arc%d_cell%d: qf, at=%e;\n'%(i_arc, i_cell, s_start_cell)
		sequence += 'sf_arc%d_cell%d: sf, at=%e;\n'%(i_arc, i_cell, s_start_cell)
		
		for i_bend in xrange(n_dip_half_cell):
			sequence += 'mb_arc%d_cell%d_%d: mb, at=%e;\n'%(i_arc, i_cell, i_bend, 
									s_start_cell+(i_bend+1)*L_halfcell/(n_dip_half_cell+1))

		sequence += 'qd_arc%d_cell%d: qd, at=%e;\n'%(i_arc, i_cell, s_start_cell+L_halfcell)
		sequence += 'sd_arc%d_cell%d: sd, at=%e;\n'%(i_arc, i_cell, s_start_cell+L_halfcell)

		for i_bend in xrange(n_dip_half_cell):
			sequence += 'mb_arc%d_cell%d_%d: mb, at=%e;\n'%(i_arc, i_cell, i_bend+n_dip_half_cell, 
									s_start_cell+(i_bend+1)*L_halfcell/(n_dip_half_cell+1)+L_halfcell)

		s_start_cell += L_halfcell*2

	#Dispesion suppressor - empty cell
	sequence += 'qf_arc%d_dsrighte: qf, at=%e;\n'%(i_arc, s_start_cell)	
	sequence += 'qd_arc%d_dsrighte: qd, at=%e;\n'%(i_arc, s_start_cell+L_halfcell)
	s_start_cell += L_halfcell*2
	#Dispesion suppressor - full cell
	sequence += 'qf_arc%d_dsrightf: qf, at=%e;\n'%(i_arc, s_start_cell)	
	for i_bend in xrange(n_dip_half_cell):
		sequence += 'mb_arc%d_dsright_%d: mb, at=%e;\n'%(i_arc, i_bend, 
								s_start_cell+(i_bend+1)*L_halfcell/(n_dip_half_cell+1))
	sequence += 'qd_arc%d_dsrightdf: qd, at=%e;\n'%(i_arc, s_start_cell+L_halfcell)
	for i_bend in xrange(n_dip_half_cell):
		sequence += 'mb_arc%d_dsright_%d: mb, at=%e;\n'%(i_arc, i_bend+n_dip_half_cell, 
								s_start_cell+(i_bend+1)*L_halfcell/(n_dip_half_cell+1)+L_halfcell)
	s_start_cell += L_halfcell*2


	#Straight
	for i_cell in xrange(n_regcells_straight/2):
		sequence += 'at_qf_ss%dL_cell%d: marker at=%e;\n'%(i_arc, i_cell, s_start_cell)
		sequence += 'qf_ss%dL_cell%d: qf, at=%e;\n'%(i_arc, i_cell, s_start_cell)
		sequence += 'qd_ss%dL_cell%d: qd, at=%e;\n'%(i_arc, i_cell, s_start_cell+L_halfcell)
		s_start_cell += L_halfcell*2

	#Insertion
	sequence += 'qf3L_ss%d: qf3, at=%e;\n'%(i_arc, s_start_cell)
	sequence += 'qd3L_ss%d: qd3, at=%e;\n'%(i_arc, s_start_cell+L_halfcell)
	s_start_cell += L_halfcell*2
	sequence += 'qf2L_ss%d: qf2, at=%e;\n'%(i_arc, s_start_cell)
	sequence += 'qd2L_ss%d: qd2, at=%e;\n'%(i_arc, s_start_cell+L_halfcell)
	s_start_cell += L_halfcell*2
	sequence += 'qf1L_ss%d: qf1, at=%e;\n'%(i_arc, s_start_cell)
	sequence += 'qd1_ss%d: qd1, at=%e;\n'%(i_arc, s_start_cell+L_halfcell)
	s_start_cell += L_halfcell*2

	sequence += 'qf1R_ss%d: qf1, at=%e;\n'%(i_arc, s_start_cell)
	sequence += 'qd2R_ss%d: qd2, at=%e;\n'%(i_arc, s_start_cell+L_halfcell)
	s_start_cell += L_halfcell*2
	sequence += 'qf2R_ss%d: qf2, at=%e;\n'%(i_arc, s_start_cell)
	sequence += 'qd3R_ss%d: qd3, at=%e;\n'%(i_arc, s_start_cell+L_halfcell)
	s_start_cell += L_halfcell*2
	sequence += 'qf3R_ss%d: qf3, at=%e;\n'%(i_arc, s_start_cell)
	sequence += 'qd4R_ss%d: qd4, at=%e;\n'%(i_arc, s_start_cell+L_halfcell)
	s_start_cell += L_halfcell*2

	for i_cell in xrange(n_regcells_straight/2):
		sequence += 'at_qf_ss%dR_cell%d: marker at=%e;\n'%(i_arc, i_cell, s_start_cell)
		sequence += 'qf_ss%d_cell%dR: qf, at=%e;\n'%(i_arc, i_cell, s_start_cell)
		sequence += 'qd_ss%d_cell%dR: qd, at=%e;\n'%(i_arc, i_cell, s_start_cell+L_halfcell)
		s_start_cell += L_halfcell*2

sequence+='''

end_machine: marker at=circum;
endsequence;

'''

with open('sequence.seq', 'w') as fid:
	fid.write(sequence)





madxscript = '''
TITLE, 'Toyring'; 
call file="sequence.seq";

!define the beam and its properties
beam, particle = proton, sequence=toyring, energy = 6500., NPART=1.05E11, sige= 2.5e;

! define the desired output
use, sequence=toyring; 


! execute the TWISS command 
!twiss,save,centre,file=twiss.out;

! match tune 
match, sequence=toyring; 
	vary,name=kqf, step=0.00001; 
	vary,name=kqd, step=0.00001; 
	global,sequence=toyring,Q1=!!Qx!!; 
	global,sequence=toyring,Q2=!!Qy!!; 
	Lmdif, calls=10, tolerance=1.0e-21;
endmatch;

! match chromaticity 
match, sequence=toyring; 
	vary,name=ksf, step=0.00001; 
	vary,name=ksd, step=0.00001; 
	global,sequence=toyring,DQ1=!!Qpx!!; 
	global,sequence=toyring,DQ2=!!Qpy!!; 
	Lmdif, calls=10, tolerance=1.0e-21;
endmatch;

! execute the TWISS command 
select,flag=twiss,column=name,s,x,y,mux,betx,
                         muy,bety,dx,dy, angle, k0l, k1l;
twiss,save,centre,file=twiss.out;

savebeta, label = leftb0, place = at_qf_ss0L_cell0;
savebeta, label = rightb0, place = at_qf_ss0R_cell0;


plot, haxis=s, vaxis=betx, bety, colour=100;
plot, haxis=s, vaxis=dx, dy, colour=100;
stop;


'''.replace('!!Qx!!', '%e'%Qx).replace('!!Qy!!', '%e'%Qy).replace('!!Qpx!!', '%e'%Qpx).replace('!!Qpy!!', '%e'%Qpy)


with open('automad.madx', 'w') as fid:
	fid.write(madxscript)

import os
os.system('../madx automad.madx')



