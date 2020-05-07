#!/usr/bin/env python
from sympy import N
from sympy import sqrt, factorial, factorial2, pi
from sympy import Symbol
l = Symbol('l', integer=True)
m = Symbol('m', integer=True)
x = Symbol('x')
######
PRECISION="double" # "single" "double" "quadruple"
LANGUAGE="Fortran" # "Fortran" "c++" "c++avx2" "c++avx512"
DO_DERIV=True
PHASE="aims" # "aims" "Condon–Shortley" None
TOTAL_LMAX=20
######

def Klm(l,m):
    return sqrt((2*l+1)*factorial(l-abs(m))/(4*pi*factorial(l+abs(m))))
def Pmm(m):
    return factorial2(2*m-1)*(-1)**m

def sConst(d):
    precision = 20
    if LANGUAGE == "Fortran":
        if PRECISION=="single":
            precision = 11
            return str(N(d,precision))+""
        if PRECISION=="double":
            precision = 20
            return str(N(d,precision))+"d0"
        if PRECISION=="quadruple":
            precision = 36
            return str(N(d,precision))+"q0"
    else:
        return str(N(d,precision))
def sMul(s1,s2):
    return s1 + "*" + s2
def sAdd(s1,s2):
    return s1 + "+" + s2
def sSub(s1,s2):
    return s1 + "-" + s2
def sAssign(sVar, sRHS):
    return sVar + " = " + sRHS
def sAssignDeriv(sVar, sRHS):
    if (DO_DERIV):
        return sVar + " = " + sRHS + ";\n"
    else:
        return ""
def sSHIndex(idx):
    if LANGUAGE=="Fortran":
        return "pSH(" + str(idx+1) + ")"
    else:
        return "pSH[" + str(idx) + "]"
def sdSHdphi_sinthetaIndex(idx):
    if LANGUAGE=="Fortran":
        return "pdSHdphi_sintheta(" + str(idx+1) + ")"
    else:
        return "pdSHdphi_sintheta[" + str(idx) + "]"
def sdSHdthetaIndex(idx):
    if LANGUAGE=="Fortran":
        return "pdSHdtheta(" + str(idx+1) + ")"
    else:
        return "pdSHdtheta[" + str(idx) + "]"

def sRuleA(m,fVal):
    return sConst((Pmm(m)*Klm(m,m)*fVal))
def sRuleB(m,fVal):
    return sMul(sConst(((2*m+1.0)*Pmm(m)*Klm(m+1,m)*fVal)),"fZ")
def sRuleC(l,m,sPm1,sPm2):
    fA=Klm(l,m)/Klm(l-1,m)*(2*l-1)/(l-m)
    fB=-Klm(l,m)/Klm(l-2,m)*(l+m-1)/(l-m)
    return sAdd(sMul(sMul(sConst(fA),"fZ"),sPm1),sMul(sConst(fB),sPm2))
def sRuleD(m,fVal):
    l = m + 2
    return sAdd(sMul(sConst(( (2*m+3)*(2*m+1)*Pmm(m)/2*Klm(l,m)*fVal )),"fZ2"),sConst((-1*(2*m+1)*  Pmm(m)/2*Klm(l,m)*fVal)))
def sRuleE(m,fVal):
    l = m + 3
    Pu = Pmm(m)
    fA = (2*m+5)*(2*m+3)*(2*m+1)*Pu/6*Klm(m+3,m)*fVal
    fB = -fVal*Klm(m+3,m)*((2*m+5)*(2*m+1)*Pu/6 + (2*m+2)*(2*m+1)*Pu/3)
    str_tmp = "(" + sAdd(sMul(sConst(fA),"fZ2"),sConst(fB)) + ")"
    return sMul("fZ",str_tmp)
def sRuleDeriv(l,m,sPml,sPml_1):
    fml = l
    fml_1 = -(l+m)*Klm(l,m)/Klm(l-1,m)
    str_ml = sMul(sMul(sConst(fml),"costheta"),sPml)
    str_ml_1 = sMul(sConst(fml_1),sPml_1)
    if m >= l:
        return str_ml
    else:
        return sAdd(str_ml,str_ml_1)
def sRuleDeriv_1(l,sPml):
    fml = Klm(l,0) / (sqrt(2) * Klm(l,1))
    return sMul(sMul(sConst(fml),"sintheta"),sPml)

def sCreateSinReccur(sCL,sSL):
    return sAdd(sMul("fX",sSL),sMul("fY",sCL))
def sCreateCosReccur(sCL,sSL):
    return sSub(sMul("fX",sCL),sMul("fY",sSL))

def comment():
    if LANGUAGE == "Fortran":
        return "!"
    else:
        return "//"
def datatype():
    if LANGUAGE == "Fortran":
        if PRECISION == "double":
            return "real*8"
        if PRECISION == "single":
            return "real*4"
        if PRECISION == "quadruple":
            return "real*16"
    else:
        return "double"

def BuildSHEvalCode(lmax):
    s_output = ""
    if LANGUAGE == "Fortran":
        if DO_DERIV:
            s_output += "subroutine SHEvalderiv" + str(lmax) + "(sintheta,costheta,sinphi,cosphi,pSH,pdSHdtheta,pdSHdphi_sintheta)\n"
        else:
            s_output += "subroutine SHEval" + str(lmax) + "(sintheta,costheta,sinphi,cosphi,pSH)\n"
        s_output += "implicit none\n"
        s_output += datatype() + ", intent(in) :: sintheta, costheta, sinphi, cosphi\n"
        s_output += datatype() + ", intent(out) :: pSH(" + str((lmax+1)*(lmax+1)) + ")\n"
        if DO_DERIV:
            s_output += datatype() + ", intent(out) :: pdSHdtheta(" + str((lmax+1)*(lmax+1)) + ")\n"
            s_output += datatype() + ", intent(out) :: pdSHdphi_sintheta(" + str((lmax+1)*(lmax+1)) + ")\n"
    else:
        if DO_DERIV:
            s_output += "void SHEval" + str(lmax) +     "(const double sintheta, const double costheta, const double sinphi, const double cosphi, double *pSH, double *pdSHdtheta, double *pdSHdphi_sintheta)\n{\n"
        else:
            s_output += "void SHEval" + str(lmax) +     "(const double sintheta, const double costheta, const double sinphi, const double cosphi, double *pSH)\n{\n"
    if (lmax > 0):
        if LANGUAGE == "Fortran":
            s_output += datatype() + " :: fX, fY, fZ\n"
            s_output += datatype() + " :: fC0,fC1,fS0,fS1,fTmpA,fTmpB,fTmpC\n"
            s_output += datatype() + " :: fC0_1,fC1_1,fS0_1,fS1_1,fTmpA_1,fTmpB_1,fTmpC_1\n"
            s_output += datatype() + " :: fTmpA_2,fTmpB_2,fTmpC_2\n"
            s_output += datatype() + " :: fZ2 \n"
        else:
            s_output += datatype() +" fX, fY, fZ;\n"
            s_output += datatype() +" fC0,fC1,fS0,fS1,fTmpA,fTmpB,fTmpC;\n"
            s_output += datatype() +" fC0_1,fC1_1,fS0_1,fS1_1,fTmpA_1,fTmpB_1,fTmpC_1;\n"
            s_output += datatype() +" fTmpA_2,fTmpB_2,fTmpC_2;\n"
        s_output += "fX = sintheta * cosphi;\n"
        s_output += "fY = sintheta * sinphi;\n"
        s_output += "fZ = costheta;\n"

    if (lmax >= 2):
        if LANGUAGE == "Fortran":
            s_output += "fZ2 = fZ*fZ;\n\n"
        else:
            s_output += datatype() + " fZ2 = fZ*fZ;\n\n"
    else:
        s_output += "\n"

    s_output += "\n"+comment()+"m = 0 \n"
    s_output += sAssign(sSHIndex(0),sConst(Klm(0,0))) + ";\n"
    s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(0),sConst(0))
    s_output += sAssignDeriv(sdSHdthetaIndex(0),sConst(0))
    if (lmax == 0):
        if LANGUAGE == "Fortran":
            s_output += "end subroutine\n\n"
        else:
            s_output += "}\n\n"
        return s_output
    m = 0
    l = 1
    idx = l*l + l
    s_output += sAssign(sSHIndex(idx), sRuleB(m,1)) + ";\n"
    s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idx),sConst(0))
    if (lmax >= 2):
        l = 2
        idx = l*l+l
        s_output += sAssign(sSHIndex(idx), sRuleD(m,1)) + ";\n"
        s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idx),sConst(0))
    if (lmax >= 3):
        l = 3
        idx = l*l+l
        s_output += sAssign(sSHIndex(idx), sRuleE(m,1)) + ";\n"
        s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idx),sConst(0))
    for l in range(4,lmax+1):
        sPm1 = sSHIndex((l-1)*(l-1)+(l-1))
        sPm2 = sSHIndex((l-2)*(l-2)+(l-2))
        idx = l*l+l
        s_output += sAssign(sSHIndex(idx), sRuleC(l,m,sPm1,sPm2)) + ";\n"
        s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idx),sConst(0))

    s_output += "\n"+comment()+ "m = 1 \n"
    if DO_DERIV:
        s_output += "fC0_1 = cosphi;\nfS0_1 = sinphi;\n"
    s_output += "fC0 = fX;\nfS0 = fY;\n"

    sC_forderiv = ["fC0_1","fC1_1"]
    sS_forderiv = ["fS0_1","fS1_1"]

    sC = ["fC0","fC1"]
    sS = ["fS0","fS1"]
    sPrev = ["fTmpA","fTmpB","fTmpC"]
    sPrev_phi = ["fTmpA_1","fTmpB_1","fTmpC_1"]
    sPrev_theta = ["fTmpA_2","fTmpB_2","fTmpC_2"]

    idxSC=0
    idxP = 0

    for m in range(1,lmax):
        l = m
        idxC = l*l+l+m
        idxS = l*l+l-m
        idxP = 0
        s_output += sPrev[idxP] + " = " + sRuleA(m,sqrt(2)) + ";\n"
        if DO_DERIV:
            s_output += sPrev_phi[idxP] + " = " + sMul(sConst(m),sPrev[idxP]) + ";\n"
            s_output += sPrev_theta[idxP] + " = " + sRuleDeriv(l,m,sPrev[idxP],sPrev[(idxP-1)%3]) + ";\n"
            if m == 1:
                s_output += sAssignDeriv(sdSHdthetaIndex(idxC-1),sRuleDeriv_1(l,sPrev[idxP]))
        if (m%2 == 1) and (PHASE=="Condon–Shortley"):
            s_output += sAssign(sSHIndex(idxC),sMul("-"+sPrev[idxP],sC[idxSC&1])) + ";\n"
        else:
            s_output += sAssign(sSHIndex(idxC),sMul(sPrev[idxP],sC[idxSC&1])) + ";\n"
        if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
            s_output += sAssign(sSHIndex(idxS),sMul("-"+sPrev[idxP],sS[idxSC&1])) + ";\n"
        else:
            s_output += sAssign(sSHIndex(idxS),sMul(sPrev[idxP],sS[idxSC&1])) + ";\n"
        if (m%2 == 1) and (PHASE=="Condon–Shortley"):
            s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxC),sMul(sPrev_phi[idxP],sS_forderiv[idxSC&1]))
        else:
            s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxC),sMul("-"+sPrev_phi[idxP],sS_forderiv[idxSC&1]))
        if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
            s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxS),sMul("-"+sPrev_phi[idxP],sC_forderiv[idxSC&1]))
        else:
            s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxS),sMul(sPrev_phi[idxP],sC_forderiv[idxSC&1]))
        if (m%2 == 1) and (PHASE=="Condon–Shortley"):
            s_output += sAssignDeriv(sdSHdthetaIndex(idxC),sMul("-"+sPrev_theta[idxP],sC_forderiv[idxSC&1]))
        else:
            s_output += sAssignDeriv(sdSHdthetaIndex(idxC),sMul(sPrev_theta[idxP],sC_forderiv[idxSC&1]))
        if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
            s_output += sAssignDeriv(sdSHdthetaIndex(idxS),sMul("-"+sPrev_theta[idxP],sS_forderiv[idxSC&1]))
        else:
            s_output += sAssignDeriv(sdSHdthetaIndex(idxS),sMul(sPrev_theta[idxP],sS_forderiv[idxSC&1]))
        if (m+1 <= lmax):
            l += 1
            idxC = l*l+l+m
            idxS = l*l+l-m
            idxP += 1
            s_output += sPrev[idxP] + " = " + sRuleB(m,sqrt(2)) + ";\n"
            if DO_DERIV:
                s_output += sPrev_phi[idxP] + " = " + sMul(sConst(m),sPrev[idxP]) + ";\n"
                s_output += sPrev_theta[idxP] + " = " + sRuleDeriv(l,m,sPrev[idxP],sPrev[(idxP-1)%3]) + ";\n"
                if m == 1:
                    s_output += sAssignDeriv(sdSHdthetaIndex(idxC-1),sRuleDeriv_1(l,sPrev[idxP]))
            if (m%2 == 1) and (PHASE=="Condon–Shortley"):
                s_output += sAssign(sSHIndex(idxC),sMul("-"+sPrev[idxP],sC[idxSC&1])) + ";\n"
            else:
                s_output += sAssign(sSHIndex(idxC),sMul(sPrev[idxP],sC[idxSC&1])) + ";\n"
            if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
                s_output += sAssign(sSHIndex(idxS),sMul("-"+sPrev[idxP],sS[idxSC&1])) + ";\n"
            else:
                s_output += sAssign(sSHIndex(idxS),sMul(sPrev[idxP],sS[idxSC&1])) + ";\n"
            if (m%2 == 1) and (PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxC),sMul(sPrev_phi[idxP],sS_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxC),sMul("-"+sPrev_phi[idxP],sS_forderiv[idxSC&1]))
            if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxS),sMul("-"+sPrev_phi[idxP],sC_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxS),sMul(sPrev_phi[idxP],sC_forderiv[idxSC&1]))
            if (m%2 == 1) and (PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdthetaIndex(idxC),sMul("-"+sPrev_theta[idxP],sC_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdthetaIndex(idxC),sMul(sPrev_theta[idxP],sC_forderiv[idxSC&1]))
            if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdthetaIndex(idxS),sMul("-"+sPrev_theta[idxP],sS_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdthetaIndex(idxS),sMul(sPrev_theta[idxP],sS_forderiv[idxSC&1]))
        if (m+2 <= lmax):
            l += 1
            idxC = l*l+l+m
            idxS = l*l+l-m
            idxP += 1
            s_output += sPrev[idxP] + " = " + sRuleD(m,sqrt(2)) + ";\n"
            if DO_DERIV:
                s_output += sPrev_phi[idxP] + " = " + sMul(sConst(m),sPrev[idxP]) + ";\n"
                s_output += sPrev_theta[idxP] + " = " + sRuleDeriv(l,m,sPrev[idxP],sPrev[(idxP-1)%3]) + ";\n"
                if m == 1:
                    s_output += sAssignDeriv(sdSHdthetaIndex(idxC-1),sRuleDeriv_1(l,sPrev[idxP]))
            if (m%2 == 1) and (PHASE=="Condon–Shortley"):
                s_output += sAssign(sSHIndex(idxC),sMul("-"+sPrev[idxP],sC[idxSC&1])) + ";\n"
            else:
                s_output += sAssign(sSHIndex(idxC),sMul(sPrev[idxP],sC[idxSC&1])) + ";\n"
            if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
                s_output += sAssign(sSHIndex(idxS),sMul("-"+sPrev[idxP],sS[idxSC&1])) + ";\n"
            else:
                s_output += sAssign(sSHIndex(idxS),sMul(sPrev[idxP],sS[idxSC&1])) + ";\n"
            if (m%2 == 1) and (PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxC),sMul(sPrev_phi[idxP],sS_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxC),sMul("-"+sPrev_phi[idxP],sS_forderiv[idxSC&1]))
            if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxS),sMul("-"+sPrev_phi[idxP],sC_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxS),sMul(sPrev_phi[idxP],sC_forderiv[idxSC&1]))
            if (m%2 == 1) and (PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdthetaIndex(idxC),"-"+sMul(sPrev_theta[idxP],sC_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdthetaIndex(idxC),sMul(sPrev_theta[idxP],sC_forderiv[idxSC&1]))
            if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdthetaIndex(idxS),sMul("-"+sPrev_theta[idxP],sS_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdthetaIndex(idxS),sMul(sPrev_theta[idxP],sS_forderiv[idxSC&1]))
        if (m+3 <= lmax):
            l += 1
            idxC = l*l+l+m
            idxS = l*l+l-m
            idxP += 1
            s_output += sPrev[idxP%3] + " = " + sRuleE(m,sqrt(2)) + ";\n"
            if DO_DERIV:
                s_output += sPrev_phi[idxP%3] + " = " + sMul(sConst(m),sPrev[idxP%3]) + ";\n"
                s_output += sPrev_theta[idxP%3] + " = " + sRuleDeriv(l,m,sPrev[idxP%3],sPrev[(idxP-1)%3]) + ";\n"
                if m == 1:
                    s_output += sAssignDeriv(sdSHdthetaIndex(idxC-1),sRuleDeriv_1(l,sPrev[idxP%3]))
            if (m%2 == 1) and (PHASE=="Condon–Shortley"):
                s_output += sAssign(sSHIndex(idxC),sMul("-"+sPrev[idxP%3],sC[idxSC&1])) + ";\n"
            else:
                s_output += sAssign(sSHIndex(idxC),sMul(sPrev[idxP%3],sC[idxSC&1])) + ";\n"
            if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
                s_output += sAssign(sSHIndex(idxS),sMul("-"+sPrev[idxP%3],sS[idxSC&1])) + ";\n"
            else:
                s_output += sAssign(sSHIndex(idxS),sMul(sPrev[idxP%3],sS[idxSC&1])) + ";\n"
            if (m%2 == 1) and (PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxC),sMul(sPrev_phi[idxP%3],sS_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxC),sMul("-"+sPrev_phi[idxP%3],sS_forderiv[idxSC&1]))
            if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxS),sMul("-"+sPrev_phi[idxP%3],sC_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxS),sMul(sPrev_phi[idxP%3],sC_forderiv[idxSC&1]))
            if (m%2 == 1) and (PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdthetaIndex(idxC),sMul("-"+sPrev_theta[idxP%3],sC_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdthetaIndex(idxC),sMul(sPrev_theta[idxP%3],sC_forderiv[idxSC&1]))
            if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdthetaIndex(idxS),sMul("-"+sPrev_theta[idxP%3],sS_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdthetaIndex(idxS),sMul(sPrev_theta[idxP%3],sS_forderiv[idxSC&1]))
        for l in range(m+4,lmax+1):
            idxC = l*l+l+m
            idxS = l*l+l-m
            idxP += 1
            s_output += sPrev[idxP%3] + " = " + sRuleC(l,m,sPrev[(idxP+3-1)%3],sPrev[(idxP+3-2)%3])   + ";\n"
            if DO_DERIV:
                s_output += sPrev_phi[idxP%3] + " = " + sMul(sConst(m),sPrev[idxP%3]) + ";\n"
                s_output += sPrev_theta[idxP%3] + " = " + sRuleDeriv(l,m,sPrev[idxP%3],sPrev[(idxP-1)%3]) + ";\n"
                if m == 1:
                    s_output += sAssignDeriv(sdSHdthetaIndex(idxC-1),sRuleDeriv_1(l,sPrev[idxP%3]))
            if (m%2 == 1) and (PHASE=="Condon–Shortley"):
                s_output += sAssign(sSHIndex(idxC),sMul("-"+sPrev[idxP%3],sC[idxSC&1])) + ";\n"
            else:
                s_output += sAssign(sSHIndex(idxC),sMul(sPrev[idxP%3],sC[idxSC&1])) + ";\n"
            if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
                s_output += sAssign(sSHIndex(idxS),sMul("-"+sPrev[idxP%3],sS[idxSC&1])) + ";\n"
            else:
                s_output += sAssign(sSHIndex(idxS),sMul(sPrev[idxP%3],sS[idxSC&1])) + ";\n"
            if (m%2 == 1) and (PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxC),sMul(sPrev_phi[idxP%3],sS_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxC),sMul("-"+sPrev_phi[idxP%3],sS_forderiv[idxSC&1]))
            if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxS),sMul("-"+sPrev_phi[idxP%3],sC_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxS),sMul(sPrev_phi[idxP%3],sC_forderiv[idxSC&1]))
            if (m%2 == 1) and (PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdthetaIndex(idxC),sMul("-"+sPrev_theta[idxP%3],sC_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdthetaIndex(idxC),sMul(sPrev_theta[idxP%3],sC_forderiv[idxSC&1]))
            if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
                s_output += sAssignDeriv(sdSHdthetaIndex(idxS),sMul("-"+sPrev_theta[idxP%3],sS_forderiv[idxSC&1]))
            else:
                s_output += sAssignDeriv(sdSHdthetaIndex(idxS),sMul(sPrev_theta[idxP%3],sS_forderiv[idxSC&1]))
        s_output += "\n"+comment()+ "m = "+str(m+1)+ "\n"
        if DO_DERIV:
            s_output += sC_forderiv[(idxSC+1)&1] + " = " + sCreateCosReccur(sC_forderiv[idxSC&1],sS_forderiv[idxSC&1]) + ";\n"
            s_output += sS_forderiv[(idxSC+1)&1] + " = " + sCreateSinReccur(sC_forderiv[idxSC&1],sS_forderiv[idxSC&1]) + ";\n"
            s_output += sC[(idxSC+1)&1] + " = " + sMul("sintheta",sC_forderiv[(idxSC+1)&1]) + ";\n"
            s_output += sS[(idxSC+1)&1] + " = " + sMul("sintheta",sS_forderiv[(idxSC+1)&1]) + ";\n"
        else:
            s_output += sC[(idxSC+1)&1] + " = " + sCreateCosReccur(sC[idxSC&1],sS[idxSC&1]) + ";\n"
            s_output += sS[(idxSC+1)&1] + " = " + sCreateSinReccur(sC[idxSC&1],sS[idxSC&1]) + ";\n"
        idxSC += 1

    m += 1
    l = lmax
    idxC = l*l+l+m
    idxS = l*l+l-m
    idxP=(idxP+1)%3

    s_output += sPrev[idxP] + " = " + sRuleA(m,sqrt(2)) + ";\n"
    if DO_DERIV:
        s_output += sPrev_phi[idxP] + " = " + sMul(sConst(m),sPrev[idxP]) + ";\n"
        s_output += sPrev_theta[idxP] + " = " + sRuleDeriv(l,m,sPrev[idxP%3],sPrev[(idxP-1)%3]) + ";\n"
        if m == 1:
            s_output += sAssignDeriv(sdSHdthetaIndex(idxC-1),sRuleDeriv_1(l,sPrev[idxP]))
    if (m%2 == 1) and (PHASE=="Condon–Shortley"):
        s_output += sAssign(sSHIndex(idxC),sMul("-"+sPrev[idxP],sC[idxSC&1])) + ";\n"
    else:
        s_output += sAssign(sSHIndex(idxC),sMul(sPrev[idxP],sC[idxSC&1])) + ";\n"
    if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
        s_output += sAssign(sSHIndex(idxS),sMul("-"+sPrev[idxP],sS[idxSC&1])) + ";\n"
    else:
        s_output += sAssign(sSHIndex(idxS),sMul(sPrev[idxP],sS[idxSC&1])) + ";\n"
    if (m%2 == 1) and (PHASE=="Condon–Shortley"):
        s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxC),sMul(sPrev_phi[idxP],sS_forderiv[idxSC&1]))
    else:
        s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxC),sMul("-"+sPrev_phi[idxP],sS_forderiv[idxSC&1]))
    if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
        s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxS),sMul("-"+sPrev_phi[idxP],sC_forderiv[idxSC&1]))
    else:
        s_output += sAssignDeriv(sdSHdphi_sinthetaIndex(idxS),sMul(sPrev_phi[idxP],sC_forderiv[idxSC&1]))
    if (m%2 == 1) and (PHASE=="Condon–Shortley"):
        s_output += sAssignDeriv(sdSHdthetaIndex(idxC),sMul("-"+sPrev_theta[idxP],sC_forderiv[idxSC&1]))
    else:
        s_output += sAssignDeriv(sdSHdthetaIndex(idxC),sMul(sPrev_theta[idxP],sC_forderiv[idxSC&1]))
    if (m%2 == 1) and (PHASE=="aims" or PHASE=="Condon–Shortley"):
        s_output += sAssignDeriv(sdSHdthetaIndex(idxS),sMul("-"+sPrev_theta[idxP],sS_forderiv[idxSC&1]))
    else:
        s_output += sAssignDeriv(sdSHdthetaIndex(idxS),sMul(sPrev_theta[idxP],sS_forderiv[idxSC&1]))
    if LANGUAGE == "Fortran":
        s_output += "end subroutine\n\n"
    else:
        s_output += "}\n\n"
    return s_output

def BuildSHEvalinterface(l_max):
    s_output_interface = ""
    if LANGUAGE == "Fortran":
        if DO_DERIV:
            s_output_interface += "subroutine SHEvalderiv(lmax,sintheta,costheta,sinphi,cosphi,pSH,pdSHdtheta,pdSHdphi_sintheta) \n"
        else:
            s_output_interface += "subroutine SHEval(lmax,sintheta,costheta,sinphi,cosphi,pSH) \n"
        s_output_interface += "implicit none \n"
        s_output_interface += "integer, intent(in) :: lmax \n"
        s_output_interface += datatype() + ", intent(in) :: sintheta,costheta,sinphi,cosphi \n"
        s_output_interface += datatype() + ", intent(inout) :: pSH((lmax+1)*(lmax+1)) \n"
        if DO_DERIV:
            s_output_interface += datatype() + ", intent(inout) :: pdSHdtheta((lmax+1)*(lmax+1)) \n"
            s_output_interface += datatype() + ", intent(inout) :: pdSHdphi_sintheta((lmax+1)*(lmax+1)) \n"
    else:
        if DO_DERIV:
            s_output_interface += "void SHEval(const int lmax, const double sintheta, const double costheta, const double sinphi, const double cosphi, double *pSH, double *pdSHdtheta, double *pdSHdphi_sintheta) {\n"
        else:
            s_output_interface += "void SHEval(const int lmax, const double sintheta, const double costheta, const double sinphi, const double cosphi, double *pSH) {\n"
    for i in range(l_max+1):
        if (i == 0):
            if LANGUAGE == "Fortran":
                s_output_interface += "if (lmax == 0)  then \n"
            else:
                s_output_interface += "if (lmax == 0)  { \n"
        else:
            if LANGUAGE == "Fortran":
                s_output_interface += "else if (lmax == " + str(i) + ")  then \n"
            else:
                s_output_interface += "} else if (lmax == " + str(i) + ")  { \n"
        if LANGUAGE == "Fortran":
            if DO_DERIV:
                s_output_interface += "call SHEvalderiv" + str(i) + "(sintheta, costheta, sinphi, cosphi, pSH, pdSHdtheta, pdSHdphi_sintheta) \n"
            else:
                s_output_interface += "call SHEval" + str(i) + "(sintheta, costheta, sinphi, cosphi, pSH) \n"
        else:
            if DO_DERIV:
                s_output_interface += "SHEval" + str(i) + "(sintheta, costheta, sinphi, cosphi, pSH, pdSHdtheta, pdSHdphi_sintheta); \n"
            else:
                s_output_interface += "SHEval" + str(i) + "(sintheta, costheta, sinphi, cosphi, pSH); \n"
        if (i==l_max):
            if LANGUAGE == "Fortran":
                s_output_interface += "end if\n"
            else:
                s_output_interface += "}\n"
    if LANGUAGE == "Fortran":
        s_output_interface += "end subroutine\n\n"
    else:
        s_output_interface += "}\n"
    return s_output_interface

def BuildSHEvalall(l_max):
    s_output_all = ""
    for i in range(l_max+1):
        s_output_all += BuildSHEvalCode(i)
    s_output_all += BuildSHEvalinterface(l_max)
    return s_output_all

print(BuildSHEvalall(TOTAL_LMAX))

