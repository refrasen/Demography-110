2+2 == 2^2 == sqrt(16)
2+2 == 3
A= diagm([1,2,3,4,5],-1)
A[1,:] = [0,1,2,3,4,5]

A

B = [1;2;3]

#space puts off to right
C = [B [1;2;3]]

C = [1 2 3]

#semicolon puts below
D = [C; [4 5 6]]
#Comma make array
E = [A, B, C]

# Access First Column of First Element
E[1][:,1]

k = [100,0,0,0,0,0]

#Multiply matrix and vector
A*k

#dot product
dot(k,k)

# Matrix Exponentiation
A^4 * k

A^3*(A*k)

#values from table 4.2 in Demography 110 book
nLx = [4770, 4726, 4712, 4698, 4681, 4662, 4637, 4604, 4561, 4503, 4421]
nfx = [0, 0, 0, .0811, .2384, .1969, .1033, .0313, .0046, .0009,0]
nBx = [0, 0, 0, 381, 1116, 918, 479, 144, 21, 4]
l0 = 1000;
n = 5;
ffab = .4877;
(nLx[1]/(2*l0))*ffab
nLx[2]/nLx[1]
function LeslieMatrix(nLx, nfx, n, l0, ffab)
    #this function generates a leslie matrix
    #n is the width of the age intervals
    #nLx is a vector containing the PYL needed to find subdiagonal
    #nfx is the vector containing ASFR needed to find first row
    N = length(nLx)-1; #this is the size of the matrix
    A = zeros((N, N)); #initializies a zero matrix
    for i in 2:(N)
        A[i, i-1] = nLx[i]/nLx[i-1];
    end
    for i in 1:(N-1)
        A[1, i] = (nLx[1]/(2*l0))*(nfx[i]+nfx[i+1]*A[i+1, i])*ffab;
    end
    A[1, N] = (nLx[1]/(2*l0))*(nfx[N]+nfx[N+1]*nLx[length(nLx)]/nLx[N])*ffab
    return A
end
B =LeslieMatrix(nLx, nfx, n, l0, ffab)
L = [9.87, 9.62, 9.39, 9.14, 8.87, 8.26, 7.45, 6.22, 4.97, 0]
F = [0, .083, .135, .198, .061, 0, 0, 0, 0, 0]

n = 10;
#A= LeslieMatrix(L, F, n, 1, .4886)
A
K = [15, 14.3, 14.5, 13.8, 12.9, 10.8, 8.5, 5.4, 2.0]
A*K
sum(K)
sum(A*K)
K = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
K5 = B*K
B^2*K
B*K5
B^3*K
B^8*K


function ASMR(D, K)
    #D is vector of deaths in age intervals
    #K is vector of Population in age intervals
    #note, if K is from a period of T = 1, K = PPYL
    #M is vector of ASMRs
    M = zeros((length(D)));
    for i in 1:length(D)
        M[i] = D[i]/K[i];
    end
    return M;
end

#naxcalc is for finding nax's when we have 0 to 1, 1 to 5, then any other equally spaced age groups
function naxcalc(M, n, inf::Bool)
    #this function finds the nax values from ASMRs and width n interval
    #M is the vector containing ASMR values (nMx)
    #n is the width of intervals after age 5
    #inf tells us if we need to caluclate inf_a_x
    nax = zeros(length(M));
    #calculates 1a0
    nax[1] =.07 + 1.7*M[1];
    nax[2] = 1.5;
    for i in 3:(length(M)-1)
        nax[i] = n/2;
    end
    #Calculating last nax value
    if inf #if it is from some age to infinity
        nax[length(M)] = 1/M[length(M)]
    else #otherwise it is another nax s
    nax[length(M)] = n/2
end
    return nax
end


function M2ql(nax, M, n, radix)
    #this function finds nqx and lx values from nax and ASMRs
    #reminder that nqx is probability of death from ages x up until x+n
    #reminder that lx is survivorship at age x
    #nax vector of nax values
    #M vector of ASMRs
    #n is a vector containing age intervals
    N = length(M); #so I don't have to keep typing "length(M)"
    nqx = zeros((N)); #initalize nqx vector, 1 good use N
    lx = zeros((N+1)); #initialize lx vector, 2 good use N
    lx[1] = radix; #we have an initial "population"
    for i in 1:N #3 good use N
        #time to calculate nqx
        nqx[i] = (n[i]*M[i])/(1+(1-nax[i])*M[i]);
        lx[i+1] = lx[i]*(1-nqx[i]);
    end
    return nqx,lx;
end

#From DEMOG 110 PS4 problem 3
DMale = [13700, 2450, 1325];
KMale = [2029000, 8282000, 10372000];
DFem = [10900, 1850, 1000];
KFem = [1943000, 7932000, 9940000];
n = [1, 4, 5]
MMale = ASMR(DMale, KMale)

MFem = ASMR(DFem, KFem)
naxMale = naxcalc(MMale, 5, false)
naxFem = naxcalc(MFem, 5, false)

nqxMale, lxMale = M2ql(naxMale, MMale, n, 1)
nqxFem, lxFem = M2ql(naxFem, MFem, n, 1)
nqxMale[1]
naxMale[1]
MMale[1]
a = 1325/10372000
b = 2.5
n = 5
n*a/(1+(1-b)*a)
nqxFem[1]
lxMale[4]
lxFem[4]
#Demonstration of how bool works
function f(b::Bool)
    return b,b
end
f(false)

1-0.4886

lxMale[4]*(1-.4886)

.5114*.9915

.4886*.9930

.4886*lxFem[4]
.50705+.48518

function BrassModel(a, b, x)
    #See table 7.3 from DEMOG 110 book
    #we will be using the spines for brass general standard
    #a and b are alpha and beta values to plug in
    #x is the age at which we want to know survivorship, should be one of the vavlues we see below in X
    #based on radix = 1
    #This will find the survivorship estimate at age x
    X = [1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80] #ages we use for brass model
    Yx = [.83659, .71963, .66572, .63295, .61005, .54619, .51316, .45500, .38287, .31502, .24966, .18157, .10741, .02120, -.08319, -.21003, -.37459, -.58184, -.86730, -1.28571] #Spines, check here for possible wrong values if you get weird values
    #using lx = 1/(1 + exp(-2a-2bYx))
    if x <=5
        i = x;
    else
        i = Int(x/5 + 4);
    end
    lx = 1/(1 + exp(-2*a-2*b*Yx[i]));
    return lx;
end

l60 = BrassModel(1.25, .93, 60)

mod(7, 5)

60/5



Yx[16]

1 - l60

l65 = BrassModel(1.25, .93, 65)

l65-l60

l60-l65

l65*5+(2.5)(l60-l65)

.8585456*5+2.5*.03326245

l70 = BrassModel(1.25, .93, 70)

l65-l70

l70*5+(2.5(l65-l70))

.80499*5+(2.5)(.05356)

l70*5+(2.5(l65-l70)) + l65*5+(2.5)(l60-l65)

4.37588+4.15884

(l60 - l70)*5 + l70*10

10*(l70+l60)/2

2.5*(l60+l65+l65+l70)

2.5*(l60+l65)

2.5*(l65+l70)

5*(l60+l70)

#The following is problem 7.6 in DEMOG110
Z55 = 0.5*log((.89658)/(1-.89658))

l80 = .47084

l55 = .89658

function Zx(lx)
    #Zx = alpha + beta*Yx
    #for Brass Model
    Zx = .5*log(lx/(1-lx));
    return Zx
end

Z55 = Zx(l55)

Z80 = Zx(l80)

.5*log(l80/(1-l80))

l80/(1-l80)

l55/(1-l55)

1-l80

1-l55

function Yxval(x)
    #finds the Brass Model Yx spine for age x
        Yx = [.83659, .71963, .66572, .63295, .61005, .54619, .51316, .45500, .38287, .31502, .24966, .18157, .10741, .02120, -.08319, -.21003, -.37459, -.58184, -.86730, -1.28571]
    if x <= 5
        return Yx[x]
    else
        return Yx[Int(x/5 + 4)]
    end
end

Yxval(1)

Yxval(60)

Yxval(55)

Yxval(80)

#A= [1 -.08319; 1 -1.28571]
Z = [1.07989; -0.05839]
alphbet = Z'/A

alphbet*A

A*alphbet

#alphbet = A\ Z # A inverse * Z

#A*alphbet

alphbet

alphbet[1] + alphbet[2]*-1.28571
function findalphbet(l,x)
    #l is a vector of lx values we have
    #x is a vector giving corresponding ages for lx
    #this finds the alpha and beta values for the BrassModel for these lx's
    Z = Zx.(l)
    A = [ones(length(l)) Yxval.(x)];
    return A\Z

end

l = [.89, .75];
x = [1, 35];
alphbet = findalphbet(l, x)
l1 = BrassModel(alphbet[1], alphbet[2], 1)
Matrix = [1 Yxval(1); 1 Yxval(35)]
Matrix*alphbet
1- l1
Zx(l[1])
Zx(l[2])



Zx.(l)

l

alphbet[1] + alphbet[2]*.83659
alphbet
BrassModel(alphbet[1], alphbet[2], 65)

BrassModel(1.25, .93, 60)

l60 = .892
