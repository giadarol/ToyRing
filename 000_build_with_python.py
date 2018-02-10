import numpy as np


L_halfcell = 50.
phase_adv_cell = np.pi/2
n_cells_arc = 24 # to have zero dispersionin the SS needs to be a multiple of 4 
n_arcs = 8
n_dip_half_cell = 3
frac_q_x = .27
frac_q_y = .295
Qpx = 15.
Qpy = 20.

Qx = n_arcs*n_cells_arc*phase_adv_cell/2./np.pi + frac_q_x
Qy = n_arcs*n_cells_arc*phase_adv_cell/2./np.pi + frac_q_y

circum = n_arcs*n_cells_arc*L_halfcell*2.

n_dip_total = n_arcs*n_cells_arc*n_dip_half_cell*2

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

circum = %e;
toyring: sequence, refer=centre, l=circum; 

'''%circum

s_start_cell = 0.
for i_arc in xrange(n_arcs):
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
                         muy,bety,dx,dy;
twiss,save,centre,file=twiss.out;
plot, haxis=s, vaxis=betx, bety, colour=100;
stop;


'''.replace('!!Qx!!', '%e'%Qx).replace('!!Qy!!', '%e'%Qy).replace('!!Qpx!!', '%e'%Qpx).replace('!!Qpy!!', '%e'%Qpy)


with open('automad.madx', 'w') as fid:
	fid.write(madxscript)

import os
os.system('../madx automad.madx')



