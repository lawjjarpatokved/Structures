m=1.0; mm=0.001*m; cm=0.01*m; km=1000.0*m; inch= 0.0254*m;ft=inch*12
KN=1.0; N=0.001*KN
Mpa= 10**3* KN/m**2; Gpa= 10**6* KN/m**2
sec=1
tonne=1*KN*(sec**2)/m; kg=0.001*tonne
lb = 4.4482216152605 * N   # 1 lbf = 4.4482216152605 N
kip = 1000.0 * lb
slug = lb * sec**2 / ft    # mass in US (F = m·a)
lbm = 0.45359237 * kg      # 1 lbm = 0.45359237 kg
psi = lb / inch**2
ksi = kip / inch**2

import libdenavit.section.database.aisc as section
import libdenavit.section.wide_flange as database
class wf_Database:
    def __init__(self,Section_name,unit=inch):

        self.section=Section_name
        self.d=section.wide_flange_database[self.section]['d']*unit
        self.tw=section.wide_flange_database[self.section]['tw']*unit
        self.bf=section.wide_flange_database[self.section]['bf']*unit
        self.tf=section.wide_flange_database[self.section]['tf']*unit
        self.A=section.wide_flange_database[self.section]['A']*(unit**2)
        self.Ix=section.wide_flange_database[self.section]['Ix']*(unit**4)
        self.Iy=section.wide_flange_database[self.section]['Iy']*(unit**4)

class WF_Database:
        def __init__(self,Section_name,unit=inch):
        
            db = database.WideFlangeDB(Section_name)
            self.d   = db.d   * unit
            self.tw  = db.tw  * unit
            self.bf  = db.bf  * unit
            self.tf  = db.tf  * unit
            self.A   = db.A   * unit**2
            self.Ix  = db.Ix  * unit**4
            self.Zx  = db.Zx  * unit**3
            self.Sx  = db.Sx  * unit**3
            self.rx  = db.rx  * unit
            self.Iy  = db.Iy  * unit**4
            self.Zy  = db.Zy  * unit**3
            self.Sy  = db.Sy  * unit**3
            self.ry  = db.ry  * unit
            self.J   = db.J   * unit**4
            self.Cw  = db.Cw  * unit**6
            self.rts = db.rts * unit
            self.ho  = db.ho  * unit
########################################################################################################################