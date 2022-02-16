#hbnet Functions for geometry, checks and filereading
import numpy
import math
#
# VECTORS AND GEOM
#
def normal(vect):

	len_vect =  math.sqrt(vect[0]*vect[0] + 
                              vect[1]*vect[1] + 
                              vect[2]*vect[2])

	norm_vect = vect/len_vect
	return norm_vect

def cross(v1,v2):
	
	b1 = (v1[1]*v2[2] - v1[2]*v2[1])
	b2 = (v1[2]*v2[0] - v1[0]*v2[2])
	b3 = (v1[0]*v2[1] - v1[1]*v2[0])

	result = numpy.array([b1,b2,b3])
	return result

def subtract(pos1,pos2):
	xx = pos1[0] - pos2[0]
        yy = pos1[1] - pos2[1]
        zz = pos1[2] - pos2[2]
	return numpy.array([xx,yy,zz])

def distance(vec):

	distA = math.sqrt(vec[0]*vec[0] + 
			  vec[1]*vec[1] + 
		          vec[2]*vec[2])
	return distA

def angle(c1,c2,c3,c4):

	v1 = subtract(c1,c2)
	v2 = subtract(c3,c2)
	v3 = subtract(c2,c3)
	v4 = subtract(c4,c3)

	cprod1,cprod2 = cross(v1,v2),cross(v3,v4)
	norm1, norm2 = normal(cprod1), normal(cprod2)
	
	dprod = numpy.dot(norm1,norm2)
  	cprod = cross(norm1,norm2)

	angle = dprod
	
	if (angle > 1.0): angle = 1.0
	if (angle < -1.0): angle = -1.0

	angle = math.acos(angle)
        angle *= (180 / math.pi)  #RADIANS to DEGREES!
	x1, y1, z1 = cprod[0], cprod[1], cprod[2]
	x2, y2, z2 = v2
	
	n = x1 * x2 + y1 * y2 + z1 * z2
	
	if (n < 0.0):
		angle = -angle
	return angle

#
# CHECK FUNCTIONS 
#
def geom_check(atm1,atm2,atm3,angles):

	vec1 = normal(subtract(atm1,atm2))
	vec2 = normal(subtract(atm3,atm2))

	ang = math.acos(numpy.dot(vec1,vec2))
	ang = ang*(180/math.pi)
	
	if ang >= angles[0] and ang <=angles[1]:
		return 1

def dist_check(atm1,atm2,lens):
	nvec = subtract(atm1,atm2)
	dist = distance(nvec)

	if dist >= lens[0] and dist <= lens[1]:
		return 1 

#def plan_check()

def ene_cal(atmH,atmN,atmO,atmC):
	#DSSP energy def
	ene = 0.084*(distance(subtract(atmN,atmO)) +
		     distance(subtract(atmC,atmH)) -
		     distance(subtract(atmO,atmH)) -
		     distance(subtract(atmC,atmN)))*332
	return ene
def elemGen(sS,bB,check):
	from hbNetClass import ArrElem
	if check == 'A':
		if sS.res == 'res':
			elem1 = elem2 = elem3 = elem4 = 0

		else:
			if sS.A1 != 0:
				elem1 = ArrElem(sS.res,sS.reN,sS.A1N,sS.AAN, sS.A1, sS.AA, sS.AH)
			else:
				elem1 = 0
			if sS.A2 != 0:
				elem2 = ArrElem(sS.res,sS.reN,sS.A2N,sS.AAN, sS.A2, sS.AA, sS.AH)
			else:
				elem2 = 0
			elem3 = elem4 = 0
		
		elem5 = ArrElem(bB.res,bB.reN,bB.AN,bB.AAN, bB.A, bB.AA, 2)

	elif check == 'D':
		if sS.res == 'res':
			elem1 = elem2 = elem3 = elem4 = 0
		elif sS.res == 'LYS':
			elem1 = ArrElem(sS.res,sS.reN,sS.D1N,sS.DDN, sS.D1, sS.DD, sS.DH)
			elem2 = ArrElem(sS.res,sS.reN,sS.D2N,sS.DDN, sS.D2, sS.DD, sS.DH)
			elem3 = ArrElem(sS.res,sS.reN,sS.AAN,sS.DDN, sS.AA, sS.DD, sS.DH)
			elem4 = 0
		elif sS.res == 'ARG':
			elem1 = ArrElem(sS.res,sS.reN,sS.D1N,sS.DDN, sS.D1, sS.DD, sS.DH)
			elem2 = ArrElem(sS.res,sS.reN,sS.D2N,sS.DDN, sS.D2, sS.DD, sS.DH)
			elem3 = ArrElem(sS.res,sS.reN,sS.A1N,sS.AAN, sS.A1, sS.AA, sS.DH)
			elem4 =	ArrElem(sS.res,sS.reN,sS.A2N,sS.AAN, sS.A2, sS.AA, sS.DH)
		else:
			if sS.D1 != 0:
                                elem1 = ArrElem(sS.res,sS.reN,sS.D1N,sS.DDN, sS.D1, sS.DD, sS.DH)
			else:
				elem1 = 0
			if sS.D2 != 0:
				elem2 = ArrElem(sS.res,sS.reN,sS.D2N,sS.DDN, sS.D2, sS.DD, sS.DH)
			else:
				elem2 = 0
			elem3 = elem4 = 0
		
		if bB.res == 'PRO':
			elem5 = 0
		else:
			elem5 = ArrElem(bB.res,bB.reN,bB.DN,bB.DDN, bB.D , bB.DD, 2)
	return [elem1,elem2,elem3,elem4,elem5]
#
# FILE READ
#
def file_read(fileName,num,tot,sslist,PEP):
	from hbNetClass import SideChain
	from hbNetClass import BackBone
	
	file = open(fileName,'r')
	
	rep = SideChain('res',0,0,0,0,0,0,0,0,0)
	bb  = BackBone('res',0,0,0,0,0)
	
	dElemList = []
	aElemList = []

	fileread = 1
	print 'File ', num, 'of ',tot
	while fileread:
                lineC = file.readline()
                tempC = lineC.split()

		if tempC[0] == 'TER':
			fileread = 0

		else:
			# record relevant structural data for residue	
                        res = tempC[3]
                        num = tempC[4]
                        atm = tempC[2]
                        
                        xyz = [ float(tempC[5]),
                                float(tempC[6]),
                                float(tempC[7])]
				
                        if res in sslist.keys():
                                ref = sslist[res]
                                rep.res = res
                                rep.reN = num
                                if atm == ref.A1:
                                        rep.A1  = xyz
                                        rep.A1N = atm

                                        if ref.A1 == ref.DD:
                                                rep.DD  = rep.A1
                                                rep.DDN = rep.A1N

                                elif atm == ref.A2:
                                        rep.A2  = xyz
                                        rep.A2N = atm
                                elif atm == ref.AA:
                                        rep.AA  = xyz
                                        rep.AAN = atm
                                elif atm == ref.D1:
                                        rep.D1  = xyz
                                        rep.D1N = atm
                                elif atm == ref.D2:
                                        rep.D2  = xyz
                                        rep.D2N = atm
                                elif atm == ref.DD:
                                        rep.DD  = xyz
                                        rep.DDN = atm

                        if atm == PEP.D:
                                bb.D  = xyz
                                bb.DN = atm
                        elif atm == PEP.DD:
                                bb.res = res
                                bb.reN = num
                                bb.DD  = xyz
                                bb.DDN = atm
                        elif atm == PEP.A:
                                # atom Oxygen last in residue 
                                # record and append
                                bb.A  = xyz
                                bb.AN = atm

                                aElems = elemGen(rep,bb,'A')
				dElems = elemGen(rep,bb,'D')
				
				for ndx in range(0,5):
					if aElems[ndx]:
						aElemList.append(aElems[ndx])
					if dElems[ndx]:
						dElemList.append(dElems[ndx])
	
				rep = SideChain('res',0,0,0,0,0,0,0,0,0)
				bb  = BackBone('res',0,0,0,0,0)

                        elif atm == PEP.AA:
                                bb.AA  = xyz
                                bb.AAN = atm
					
	Darray = numpy.array(dElemList)
	Aarray = numpy.array(aElemList)	
	file.close()
	return numpy.transpose(Darray),numpy.transpose(Aarray)
