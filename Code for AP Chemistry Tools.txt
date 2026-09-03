from math import *

# =========================
# AP CHEM TOOLKIT V3
# TI-84 EVO
# TI-84 PLUS CE PYTHON
# ANY TI THAT HAS PYTHON CAPABILITIES
# =========================

R=0.082057
RJ=8.314
NA=6.022e23
F=96485
KW=1.0e-14

# Atomic masses (g/mol)
AM={
"H":1.008,"He":4.0026,
"Li":6.94,"Be":9.0122,"B":10.81,"C":12.011,
"N":14.007,"O":15.999,"F":18.998,"Ne":20.180,
"Na":22.990,"Mg":24.305,"Al":26.982,"Si":28.085,
"P":30.974,"S":32.06,"Cl":35.45,"Ar":39.948,
"K":39.098,"Ca":40.078,"Sc":44.956,"Ti":47.867,
"V":50.942,"Cr":52.00,"Mn":54.938,"Fe":55.845,
"Co":58.933,"Ni":58.693,"Cu":63.546,"Zn":65.38,
"Ga":69.723,"Ge":72.630,"As":74.922,"Se":78.971,
"Br":79.904,"Kr":83.798,"Rb":85.468,"Sr":87.62,
"Y":88.906,"Zr":91.224,"Nb":92.906,"Mo":95.95,
"Tc":98,"Ru":101.07,"Rh":102.91,"Pd":106.42,
"Ag":107.87,"Cd":112.41,"In":114.82,"Sn":118.71,
"Sb":121.76,"Te":127.60,"I":126.90,"Xe":131.29,
"Cs":132.91,"Ba":137.33,"La":138.91,"Ce":140.12,
"Pr":140.91,"Nd":144.24,"Pm":145,"Sm":150.36,
"Eu":151.96,"Gd":157.25,"Tb":158.93,"Dy":162.50,
"Ho":164.93,"Er":167.26,"Tm":168.93,"Yb":173.05,
"Lu":174.97,"Hf":178.49,"Ta":180.95,"W":183.84,
"Re":186.21,"Os":190.23,"Ir":192.22,"Pt":195.08,
"Au":196.97,"Hg":200.59,"Tl":204.38,"Pb":207.2,
"Bi":208.98,"Po":209,"At":210,"Rn":222,
"Fr":223,"Ra":226,"Ac":227,"Th":232.04,
"Pa":231.04,"U":238.03,"Np":237,"Pu":244,
"Am":243,"Cm":247,"Bk":247,"Cf":251,
"Es":252,"Fm":257,"Md":258,"No":259,
"Lr":266,"Rf":267,"Db":268,"Sg":269,
"Bh":270,"Hs":277,"Mt":278,"Ds":281,
"Rg":282,"Cn":285,"Nh":286,"Fl":289,
"Mc":290,"Lv":293,"Ts":294,"Og":294
}

# =========================
# BASIC HELPERS
# =========================

def pause():
    input("ENTER")

def add_comp(a,b,mult):
    for e in b:
        if e in a:
            a[e]+=b[e]*mult
        else:
            a[e]=b[e]*mult

def get_num(s,i):
    num=""
    while i<len(s) and s[i]>="0" and s[i]<="9":
        num+=s[i]
        i+=1

    if num=="":
        return 1,i

    return int(num),i

# =========================
# CHEMICAL FORMULA PARSER
# =========================

def parse_group(s,i=0):
    comp={}

    while i<len(s):

        ch=s[i]

        if ch==")":
            return comp,i+1

        elif ch=="(":
            sub,i=parse_group(s,i+1)
            mult,i=get_num(s,i)
            add_comp(comp,sub,mult)

        elif ch>="A" and ch<="Z":
            symbol=ch
            i+=1

            if i<len(s):
                if s[i]>="a" and s[i]<="z":
                    symbol+=s[i]
                    i+=1

            if symbol not in AM:
                raise ValueError("Unknown element")

            mult,i=get_num(s,i)

            if symbol in comp:
                comp[symbol]+=mult
            else:
                comp[symbol]=mult

        else:
            raise ValueError("Bad formula")

    return comp,i

def parse_formula(formula):
    formula=formula.replace(" ","")
    total={}

    # Hydrates:
    # CuSO4.5H2O
    parts=formula.split(".")

    for part in parts:

        if part=="":
            raise ValueError("Bad formula")

        i=0
        coeff=""

        while i<len(part) and part[i]>="0" and part[i]<="9":
            coeff+=part[i]
            i+=1

        if coeff=="":
            mult=1
        else:
            mult=int(coeff)

        sub=part[i:]

        if sub=="":
            raise ValueError("Bad formula")

        comp,end=parse_group(sub,0)

        if end!=len(sub):
            raise ValueError("Bad formula")

        add_comp(total,comp,mult)

    return total

def formula_mass(formula):
    comp=parse_formula(formula)
    mass=0

    for e in comp:
        mass+=AM[e]*comp[e]

    return mass

# =========================
# MOLAR MASS
# =========================

def molarmass():
    print("\nMOLAR MASS")
    f=input("Formula: ")

    try:
        m=formula_mass(f)
        print("MM =",round(m,4),"g/mol")
    except:
        print("Invalid formula")

    pause()

# =========================
# GRAMS / MOLES
# =========================

def gramsmoles():
    print("\nGRAMS / MOLES")
    f=input("Formula: ")

    try:
        mm=formula_mass(f)

        print("MM =",round(mm,4))
        print("1 g -> mol")
        print("2 mol -> g")

        c=input("> ")

        if c=="1":
            g=float(input("Grams: "))
            print("mol =",g/mm)

        elif c=="2":
            mol=float(input("Moles: "))
            print("g =",mol*mm)

    except:
        print("Invalid input")

    pause()

# =========================
# PERCENT COMPOSITION
# =========================

def percentcomp():
    print("\nPERCENT COMPOSITION")
    f=input("Formula: ")

    try:
        comp=parse_formula(f)
        mm=formula_mass(f)

        print("MM =",round(mm,4))

        for e in comp:
            mass=AM[e]*comp[e]
            pct=100*mass/mm
            print(e,round(pct,3),"%")

    except:
        print("Invalid formula")

    pause()

# =========================
# EMPIRICAL FORMULA
# =========================

def empirical():
    print("\nEMPIRICAL FORMULA")
    print("Enter masses or %")
    n=int(input("# elements: "))

    elems=[]
    moles=[]

    for i in range(n):
        e=input("Element: ")

        if e not in AM:
            print("Unknown element")
            pause()
            return

        amt=float(input("Mass/%: "))
        elems.append(e)
        moles.append(amt/AM[e])

    smallest=min(moles)

    ratios=[]
    for x in moles:
        ratios.append(x/smallest)

    # Find multiplier that makes
    # ratios closest to integers
    best_mult=1
    best_err=999

    for mult in range(1,13):
        err=0

        for x in ratios:
            y=x*mult
            err+=abs(y-round(y))

        if err<best_err:
            best_err=err
            best_mult=mult

    nums=[]

    for x in ratios:
        nums.append(int(round(x*best_mult)))

    print("Empirical formula:")

    formula=""

    for i in range(n):
        formula+=elems[i]

        if nums[i]!=1:
            formula+=str(nums[i])

    print(formula)

    print("Ratios:")
    for i in range(n):
        print(elems[i],nums[i])

    pause()

# =========================
# MOLECULAR FORMULA
# =========================

def molecular():
    print("\nMOLECULAR FORMULA")

    ef=input("Empirical formula: ")

    try:
        emm=formula_mass(ef)
        realmm=float(input("Molar mass: "))

        mult=int(round(realmm/emm))

        if mult<1:
            print("Invalid molar mass")
            pause()
            return

        comp=parse_formula(ef)

        formula=""

        for e in comp:
            num=comp[e]*mult
            formula+=e

            if num!=1:
                formula+=str(num)

        print("EF mass =",round(emm,4))
        print("Multiplier =",mult)
        print("Molecular:")
        print(formula)

    except:
        print("Invalid input")

    pause()

# =========================
# PARTICLES
# =========================

def particles():
    print("\nMOLES / PARTICLES")
    print("1 mol -> particles")
    print("2 particles -> mol")

    c=input("> ")

    try:
        if c=="1":
            mol=float(input("Moles: "))
            print("Particles =",mol*NA)

        elif c=="2":
            p=float(input("Particles: "))
            print("Moles =",p/NA)

    except:
        print("Invalid input")

    pause()

# =========================
# STOICHIOMETRY
# =========================

def stoich():
    print("\nSTOICHIOMETRY")
    print("A -> B")

    fa=input("Formula A: ")
    fb=input("Formula B: ")

    try:
        mma=formula_mass(fa)
        mmb=formula_mass(fb)

        ca=float(input("Coeff A: "))
        cb=float(input("Coeff B: "))
        ga=float(input("Grams A: "))

        mola=ga/mma
        molb=mola*cb/ca
        gb=molb*mmb

        print("mol A =",mola)
        print("mol B =",molb)
        print("g B =",gb)

    except:
        print("Invalid input")

    pause()

# =========================
# LIMITING REACTANT
# =========================

def limiting():
    print("\nLIMITING REACTANT")

    fa=input("Reactant A: ")
    fb=input("Reactant B: ")
    fp=input("Product: ")

    try:
        mma=formula_mass(fa)
        mmb=formula_mass(fb)
        mmp=formula_mass(fp)

        ca=float(input("Coeff A: "))
        cb=float(input("Coeff B: "))
        cp=float(input("Coeff P: "))

        ga=float(input("Grams A: "))
        gb=float(input("Grams B: "))

        mola=ga/mma
        molb=gb/mmb

        ra=mola/ca
        rb=molb/cb

        if ra<rb:
            print("Limiting:",fa)
            extent=ra

            excessmol=molb-extent*cb
            print("Excess",fb)
            print(excessmol*mmb,"g")

        elif rb<ra:
            print("Limiting:",fb)
            extent=rb

            excessmol=mola-extent*ca
            print("Excess",fa)
            print(excessmol*mma,"g")

        else:
            print("Perfect ratio")
            extent=ra

        pmol=extent*cp
        pg=pmol*mmp

        print("Theo yield:")
        print(pg,"g",fp)

    except:
        print("Invalid input")

    pause()

# =========================
# PERCENT YIELD
# =========================

def percentyield():
    print("\nPERCENT YIELD")

    try:
        actual=float(input("Actual g: "))
        theo=float(input("Theo g: "))

        print("% yield =",
              100*actual/theo)

    except:
        print("Invalid input")

    pause()

# =========================
# MOLARITY
# =========================

def molarity():
    print("\nMOLARITY")

    try:
        mol=float(input("Moles: "))
        L=float(input("Liters: "))

        print("M =",mol/L)

    except:
        print("Invalid input")

    pause()

# =========================
# DILUTION
# =========================

def dilution():
    print("\nM1V1 = M2V2")
    print("0 = unknown")

    try:
        m1=float(input("M1: "))
        v1=float(input("V1: "))
        m2=float(input("M2: "))
        v2=float(input("V2: "))

        zeros=0
        if m1==0: zeros+=1
        if v1==0: zeros+=1
        if m2==0: zeros+=1
        if v2==0: zeros+=1

        if zeros!=1:
            print("Enter ONE unknown")

        elif m1==0:
            print("M1 =",m2*v2/v1)

        elif v1==0:
            print("V1 =",m2*v2/m1)

        elif m2==0:
            print("M2 =",m1*v1/v2)

        elif v2==0:
            print("V2 =",m1*v1/m2)

    except:
        print("Invalid input")

    pause()

# =========================
# IDEAL GAS
# =========================

def idealgas():
    print("\nPV = nRT")
    print("atm, L, mol, K")
    print("0 = unknown")

    try:
        p=float(input("P: "))
        v=float(input("V: "))
        n=float(input("n: "))
        t=float(input("T: "))

        zeros=0
        if p==0: zeros+=1
        if v==0: zeros+=1
        if n==0: zeros+=1
        if t==0: zeros+=1

        if zeros!=1:
            print("Enter ONE unknown")

        elif p==0:
            print("P =",n*R*t/v,"atm")

        elif v==0:
            print("V =",n*R*t/p,"L")

        elif n==0:
            print("n =",p*v/(R*t),"mol")

        elif t==0:
            print("T =",p*v/(n*R),"K")

    except:
        print("Invalid input")

    pause()

# =========================
# COMBINED GAS
# =========================

def combinedgas():
    print("\nP1V1/T1=P2V2/T2")
    print("0 = unknown")

    try:
        p1=float(input("P1: "))
        v1=float(input("V1: "))
        t1=float(input("T1 K: "))
        p2=float(input("P2: "))
        v2=float(input("V2: "))
        t2=float(input("T2 K: "))

        if p1==0:
            print("P1 =",p2*v2*t1/(v1*t2))

        elif v1==0:
            print("V1 =",p2*v2*t1/(p1*t2))

        elif t1==0:
            print("T1 =",p1*v1*t2/(p2*v2))

        elif p2==0:
            print("P2 =",p1*v1*t2/(v2*t1))

        elif v2==0:
            print("V2 =",p1*v1*t2/(p2*t1))

        elif t2==0:
            print("T2 =",p2*v2*t1/(p1*v1))

    except:
        print("Invalid input")

    pause()

# =========================
# ACID / BASE
# =========================

def acidbase():
    while True:

        print("\nACID / BASE")
        print("1 pH -> H+")
        print("2 H+ -> pH")
        print("3 pOH -> OH-")
        print("4 OH- -> pOH")
        print("5 pH <-> pOH")
        print("6 Ka/Kb tools")
        print("7 Buffer")
        print("8 Back")

        c=input("> ")

        try:
            if c=="1":
                ph=float(input("pH: "))
                print("[H+] =",10**(-ph))
                pause()

            elif c=="2":
                h=float(input("[H+]: "))
                print("pH =",-log10(h))
                pause()

            elif c=="3":
                poh=float(input("pOH: "))
                print("[OH-] =",10**(-poh))
                pause()

            elif c=="4":
                oh=float(input("[OH-]: "))
                print("pOH =",-log10(oh))
                pause()

            elif c=="5":
                print("1 pH -> pOH")
                print("2 pOH -> pH")
                x=input("> ")

                if x=="1":
                    ph=float(input("pH: "))
                    print("pOH =",14-ph)

                elif x=="2":
                    poh=float(input("pOH: "))
                    print("pH =",14-poh)

                pause()

            elif c=="6":
                kakb()

            elif c=="7":
                buffer()

            elif c=="8":
                break

        except:
            print("Invalid input")
            pause()

# =========================
# Ka / Kb
# =========================

def kakb():
    while True:

        print("\nKa / Kb")
        print("1 Ka -> pKa")
        print("2 pKa -> Ka")
        print("3 Kb -> pKb")
        print("4 pKb -> Kb")
        print("5 Ka -> Kb")
        print("6 Kb -> Ka")
        print("7 Weak acid pH")
        print("8 Weak base pH")
        print("9 Back")

        c=input("> ")

        try:
            if c=="1":
                ka=float(input("Ka: "))
                print("pKa =",-log10(ka))
                pause()

            elif c=="2":
                pka=float(input("pKa: "))
                print("Ka =",10**(-pka))
                pause()

            elif c=="3":
                kb=float(input("Kb: "))
                print("pKb =",-log10(kb))
                pause()

            elif c=="4":
                pkb=float(input("pKb: "))
                print("Kb =",10**(-pkb))
                pause()

            elif c=="5":
                ka=float(input("Ka: "))
                print("Kb =",KW/ka)
                pause()

            elif c=="6":
                kb=float(input("Kb: "))
                print("Ka =",KW/kb)
                pause()

            elif c=="7":
                ka=float(input("Ka: "))
                C=float(input("Initial M: "))

                # x^2/(C-x)=Ka
                x=(-ka+sqrt(ka*ka+4*ka*C))/2

                print("[H+] =",x)
                print("pH =",-log10(x))
                pause()

            elif c=="8":
                kb=float(input("Kb: "))
                C=float(input("Initial M: "))

                x=(-kb+sqrt(kb*kb+4*kb*C))/2

                poh=-log10(x)

                print("[OH-] =",x)
                print("pOH =",poh)
                print("pH =",14-poh)
                pause()

            elif c=="9":
                break

        except:
            print("Invalid input")
            pause()

# =========================
# BUFFER
# =========================

def buffer():
    print("\nBUFFER")
    print("Henderson-Hasselbalch")
    print("pH=pKa+log(B/A)")
    print("1 Find pH")
    print("2 Find base/acid ratio")

    c=input("> ")

    try:
        if c=="1":
            pka=float(input("pKa: "))
            base=float(input("[Base]: "))
            acid=float(input("[Acid]: "))

            ph=pka+log10(base/acid)

            print("pH =",ph)

        elif c=="2":
            ph=float(input("pH: "))
            pka=float(input("pKa: "))

            ratio=10**(ph-pka)

            print("Base/Acid =",ratio)

    except:
        print("Invalid input")

    pause()

# =========================
# EQUILIBRIUM / Q
# =========================

def equilibrium():
    print("\nREACTION QUOTIENT")
    print("Q = products/reactants")
    print("Do not enter solids/liquids")

    try:
        np=int(input("# products: "))
        top=1

        for i in range(np):
            conc=float(input("Prod conc: "))
            power=float(input("Coeff: "))
            top*=conc**power

        nr=int(input("# reactants: "))
        bottom=1

        for i in range(nr):
            conc=float(input("React conc: "))
            power=float(input("Coeff: "))
            bottom*=conc**power

        q=top/bottom

        print("Q =",q)

        k=float(input("K (0 skip): "))

        if k>0:
            if q<k:
                print("Shifts RIGHT")
            elif q>k:
                print("Shifts LEFT")
            else:
                print("At equilibrium")

    except:
        print("Invalid input")

    pause()

# =========================
# CALORIMETRY
# =========================

def qcalc():
    print("\nq = mcDT")
    print("0 = unknown")

    try:
        q=float(input("q J: "))
        m=float(input("mass g: "))
        c=float(input("c J/gC: "))
        dt=float(input("Delta T: "))

        if q==0:
            print("q =",m*c*dt,"J")

        elif m==0:
            print("m =",q/(c*dt),"g")

        elif c==0:
            print("c =",q/(m*dt))

        elif dt==0:
            print("Delta T =",q/(m*c))

    except:
        print("Invalid input")

    pause()

# =========================
# GIBBS FREE ENERGY
# =========================

def gibbs():
    print("\nDG = DH - TDS")

    try:
        dh=float(input("DH kJ/mol: "))
        t=float(input("T K: "))
        ds=float(input("DS J/molK: "))

        dg=dh-t*(ds/1000)

        print("DG =",dg,"kJ/mol")

        if dg<0:
            print("Spontaneous")
        elif dg>0:
            print("Nonspontaneous")
        else:
            print("Equilibrium")

    except:
        print("Invalid input")

    pause()

# =========================
# BEER-LAMBERT
# =========================

def beer():
    print("\nBEER-LAMBERT")
    print("A = ebc")
    print("0 = unknown")

    try:
        A=float(input("A: "))
        e=float(input("e: "))
        b=float(input("b cm: "))
        c=float(input("c M: "))

        if A==0:
            print("A =",e*b*c)

        elif e==0:
            print("e =",A/(b*c))

        elif b==0:
            print("b =",A/(e*c),"cm")

        elif c==0:
            print("c =",A/(e*b),"M")

    except:
        print("Invalid input")

    pause()

# =========================
# KINETICS
# =========================

def kinetics():
    while True:

        print("\nKINETICS")
        print("1 Zero order")
        print("2 First order")
        print("3 Second order")
        print("4 Half-life")
        print("5 Back")

        c=input("> ")

        try:
            if c=="1":
                print("[A]=[A]0-kt")
                a0=float(input("[A]0: "))
                k=float(input("k: "))
                t=float(input("time: "))

                print("[A] =",a0-k*t)
                pause()

            elif c=="2":
                print("[A]=[A]0 e^-kt")
                a0=float(input("[A]0: "))
                k=float(input("k: "))
                t=float(input("time: "))

                print("[A] =",
                      a0*exp(-k*t))
                pause()

            elif c=="3":
                print("1/[A]=1/[A]0+kt")
                a0=float(input("[A]0: "))
                k=float(input("k: "))
                t=float(input("time: "))

                a=1/(1/a0+k*t)

                print("[A] =",a)
                pause()

            elif c=="4":
                print("1 Zero")
                print("2 First")
                print("3 Second")

                order=input("> ")

                a0=float(input("[A]0: "))
                k=float(input("k: "))

                if order=="1":
                    print("t1/2 =",a0/(2*k))

                elif order=="2":
                    print("t1/2 =",log(2)/k)

                elif order=="3":
                    print("t1/2 =",1/(k*a0))

                pause()

            elif c=="5":
                break

        except:
            print("Invalid input")
            pause()

# =========================
# ELECTROCHEMISTRY
# =========================

def electro():
    while True:

        print("\nELECTROCHEM")
        print("1 Nernst")
        print("2 DG = -nFE")
        print("3 E from DG")
        print("4 Eo from K")
        print("5 K from Eo")
        print("6 Back")

        c=input("> ")

        try:
            if c=="1":
                eo=float(input("Eo V: "))
                t=float(input("T K: "))
                n=float(input("n e-: "))
                q=float(input("Q: "))

                e=eo-(RJ*t/(n*F))*log(q)

                print("E =",e,"V")
                pause()

            elif c=="2":
                n=float(input("n e-: "))
                e=float(input("E V: "))

                dg=-n*F*e

                print("DG =",dg,"J/mol")
                print("DG =",dg/1000,"kJ/mol")
                pause()

            elif c=="3":
                dg=float(input("DG kJ/mol: "))
                n=float(input("n e-: "))

                e=-(dg*1000)/(n*F)

                print("E =",e,"V")
                pause()

            elif c=="4":
                k=float(input("K: "))
                n=float(input("n e-: "))
                t=float(input("T K: "))

                eo=(RJ*t/(n*F))*log(k)

                print("Eo =",eo,"V")
                pause()

            elif c=="5":
                eo=float(input("Eo V: "))
                n=float(input("n e-: "))
                t=float(input("T K: "))

                k=exp(n*F*eo/(RJ*t))

                print("K =",k)
                pause()

            elif c=="6":
                break

        except:
            print("Invalid input")
            pause()

# =========================
# CONSTANTS
# =========================

def constants():
    print("\nCHEM CONSTANTS")
    print("R=.082057 L atm/molK")
    print("R=8.314 J/molK")
    print("NA=6.022e23")
    print("F=96485 C/mol e-")
    print("Kw=1.0e-14 at 25C")
    print("STP=273.15 K, 1 atm")

    pause()

# =========================
# FORMULA MENU
# =========================

def formula_menu():
    while True:

        print("\nFORMULA TOOLS")
        print("1 Molar mass")
        print("2 Grams/Moles")
        print("3 Percent comp")
        print("4 Empirical")
        print("5 Molecular")
        print("6 Stoichiometry")
        print("7 Limiting")
        print("8 Percent yield")
        print("9 Back")

        c=input("> ")

        if c=="1":
            molarmass()
        elif c=="2":
            gramsmoles()
        elif c=="3":
            percentcomp()
        elif c=="4":
            empirical()
        elif c=="5":
            molecular()
        elif c=="6":
            stoich()
        elif c=="7":
            limiting()
        elif c=="8":
            percentyield()
        elif c=="9":
            break

# =========================
# SOLUTION / GAS MENU
# =========================

def solution_menu():
    while True:

        print("\nSOLUTIONS / GAS")
        print("1 Molarity")
        print("2 Dilution")
        print("3 Ideal gas")
        print("4 Combined gas")
        print("5 Mol/Particles")
        print("6 Back")

        c=input("> ")

        if c=="1":
            molarity()
        elif c=="2":
            dilution()
        elif c=="3":
            idealgas()
        elif c=="4":
            combinedgas()
        elif c=="5":
            particles()
        elif c=="6":
            break

# =========================
# MAIN MENU
# =========================

while True:

    print("\n=== AP CHEM V3 ===")
    print("1 Formula Tools")
    print("2 Solutions/Gas")
    print("3 Acid/Base")
    print("4 Equilibrium/Q")
    print("5 Calorimetry")
    print("6 Gibbs Energy")
    print("7 Beer-Lambert")
    print("8 Kinetics")
    print("9 Electrochem")
    print("A Constants")
    print("Q Quit")

    choice=input("> ")

    if choice=="1":
        formula_menu()

    elif choice=="2":
        solution_menu()

    elif choice=="3":
        acidbase()

    elif choice=="4":
        equilibrium()

    elif choice=="5":
        qcalc()

    elif choice=="6":
        gibbs()

    elif choice=="7":
        beer()

    elif choice=="8":
        kinetics()

    elif choice=="9":
        electro()

    elif choice=="A" or choice=="a":
        constants()

    elif choice=="Q" or choice=="q":
        break