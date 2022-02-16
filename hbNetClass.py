class ArrElem:
	def __init__(self,res,num,atm1,atm2,pos1,pos2,hyb):
		self.res   = res
		self.reN   = num
		self.nam1  = atm1
		self.nam2  = atm2
		self.crd1  = pos1
		self.crd2  = pos2
		self.ahyb  = hyb

class BackBone:
	def __init__(self,res,num,a,aa,d,dd):
		self.res = res
		self.reN = num
                self.A   = a
                self.AA  = aa
                self.D   = d
                self.DD  = dd
		
class SideChain:
	def __init__(self,res,num,a1,a2,aa,d1,d2,dd,ahyb,dhyb):
		self.res = res
		self.reN = num
		self.A1  = a1
		self.A2  = a2
		self.AA  = aa
		self.D1  = d1
		self.D2  = d2
		self.DD  = dd
		self.AH  = ahyb
		self.DH  = dhyb

