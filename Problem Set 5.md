Renee Senining

Demography 110

# Problem Set 5

#### 3.8 

B = 100,000 /year, starting age 75 till death, $e_{max} = e_{100} = 2.1$, $i = .05$ 

Using: $_nL_x = \frac{n}{2} \times (l_x+ l_{x+n})$ 

and $P = B \times (\frac{_{15}L_{75}}{1.05^{75 - 25 + 15/2}} + \frac{_{10}L_{90}}{1.05^{90 - 25 + 10/2}} + \frac{_{\infin}L_{100}}{1.05^{100 - 25 + 2.1}} )$ = $49,097.77

then, since we have to take into consideration that not everyone lives to the age 25, we multiply P by $l_{25} = {_{25}}P_0$, we end up with $47,905.18 (noting that our radix was 1, so we didn't need to do much more)



#### 8.1 

$_5q_0 = 0.185$ with $s^{HIV} = .01$, $s^{diarrhoea} = .12$, $s^{malaria} = .23$, $s^{pneumonia} = .17$, $s^{injuries} = .03$, $s^{prematurity} = .1$, $s ^{neonatal} = .17$, and $s^{residual} = .16$ 

$_5q_0^{malaria} = s^{malaria}\times{_5}q_0$ = 0.042550

$_5q_0^{other^*} = 1 - (1-{_5}q_0)^{1 - s^{malaria}}$  = 0.14573



#### 8.2 

$_5q_{90}^{A^*} = 0.097$ and $_5q_{90}^{B^*} = 0.132$ , A is cancer, B is all other causes.

$_5q_{90} = 0.097 + 0.132 - (0.097)(0.132) = .21620  $ 

$_nq_x^{s^*} = 1 - (1- {_nq_x})^{s} $

$(1-{_n}q_x)^s = 1 - {_n}q_x^{s^*} $

$s\times log(1 - {_n}q_x) = log(1 - {_nq_x^{s^*}}) $

$ s = \frac{log(1-{_n}q_x^{s^*})}{log(1-{_n}q_x)}$

$s^A = .41885 $, $s^B = .58113$ and note, when you add them up, it's very close to 1, since B = the complement to A in terms of causes of death.



#### Problem #3 

$l_0 = 1$, $\alpha = 0, \beta = 1$ 

Then I plugged these into the equation underneath table 7.3 with the corresponding spines 

$l_{20} = 0.71300 $

$l_{40} = 0.58980 $

$l_{60} = 0.39650 $



#### 9.1 

The table only gives me values from 15 and on, so I'm assuming that for age groups 0-4, 5-9, 10-14, PEM = 0.

So I will calculate the sum of each $5\times{_5}S_x$ for ages 15 to whatever age the $PEM_{ult}$ is and add that to 15 to get SMAFM.

SMAFM = 21.38165

My $PEM_{ult} = .959$, so I only used data from x = 15 to x = 45

to get $_nS_x$ values:

| 15     | 20     | 25     | 30     | 35     | 40     | 45   |
| ------ | ------ | ------ | ------ | ------ | ------ | ---- |
| .80292 | .34828 | .09072 | .02920 | .00417 | .00104 | 0.0  |

#### 10.3

$r = 0.027205$ 

$ \frac{_{10}K_{15}}{_{35}K_{15}} =  $ $\frac{_{10}L_{15}\times e^{-r \times 15}}{_{35}L_{15}\times e^{-r \times 15}}$ 

because B, $l_0$ are constants, and we can even cancel out $e^{-r\times 15}$ since r is constant.

and with $_{10}L_{15} = 659345$ and $_{35}L_{15} = 2008629$ 

**= 0.32826**

#### 10.4

 $b = 0.050360 $, find proportion aged 15 to 50 among all women in the stable population.

$\frac{_{35}K_{15}}{_{\infin}K_0}$ = $b \times \frac{_{35}L_{15}}{l_0}\times e^{-r \times 15} $

 **= 0.67261**

#### 10.5 

| Country   | K(2012) | CBR  | **$e_0$** | NRR(2012) |
| --------- | ------- | ---- | --------- | --------- |
| Indonesia | 245 mil | .019 | 71        | 1.1       |
| Pakistan  | 188 mil | .028 | 63        | 1.39      |
| Nigeria   | 170 mil | .040 | 47        | 2.21      |
| Brazil    | 194 mil | .016 | 73        | .86       |

$K_{ult}^{Indonesia}$ $\approx$ $K_{2012}^{Indonesia} \times \frac{b_{2012}\times e_0}{\sqrt{NRR}}$ 

Doing this for the other countries as well we get:

$K_{ult}^{Indonesia} = 315.12 $ million

$K_{ult}^{Pakistan} = 281.29 $ million

$K_{ult}^{Nigeria} = 215.99 $ million

The idea is that birth rates fall and we go from a growing stable population to a stationary one in a Keyfitz scenario, but Brazil's NRR < 1, which implies r < 0 in the long run. We would have birth rates increase for r = 0. (this comes from the first few paragraphs in 10.8)



The formula relies on the expectation that, because log of births should be ***falling*** by log(NRR), log(B(t)) after the drop to be bigger than log(B(+$\epsilon$ ) but to not rise back as high as the maximum log(B(-$\epsilon$)) under the old rate. The approximation supposes that log(B(t)) settles down halfway between log(B(+$\epsilon$)) and log(B(-$\epsilon$)), but with Brazil, there is no drop, the log(B(t)) with t being after the change in rates would need to be higher than any log(B(t)) before the change, since we'd need birth rates to increase for the NRR = 1.

#### Problem #6

We know that $_nK_{x+n}$ after n years would = $_nK_x\times\frac{_{n}L_{x+n}}{_{n}L_{x}}$ 

and knowing $_np_x$ = $\frac{l_{x+n}}{l_x}$ and using the approximation $_nL_x = $ $\frac{n}{2}\times(l_x+l_{x+n})$ 

$\frac{_nL_{x+n}}{_nL_x} = \frac{l_{x+n}+l_{x+2n}}{l_{x} + l_{x+n}}$ = $ \frac{l_{x+n}+l_{x+2n}}{l_{x} + l_{x+n}} \times 1 =   \frac{l_{x+n}+l_{x+2n}}{l_{x} + l_{x+n}} \times \frac{\frac{1}{l_{x+2n}}}{\frac{1}{l_{x+2n}}}$ $ = \frac{_np_{x+n}+1}{_{2n}p_x + _np_{x+n}}$

and survival probabilities multiply so we have, with n = 5, $_{10}p_x = {_5}p_x\times{_5}p_{x+5}$ 

and we are given all survival probabilities $_5p_x$ with $x \ge10$ , so we can get any $_{10}p_x$ for any $x \ge10$ 

so we can get, after 5 years, the population $_5K_{15}^{5 years later}$ = $_5K_{10} \times \frac{_5p_{15}+1}{_{10}p_{10}+{_5}p_{15}}$ 

after 10 years, the population $_5K_{20}^{10 years later}$ = $_5K_{15}^{5 years later} \times \frac{_5p_{20}+1}{_{10}p_{15}+{_5}p_{20}}$ 

and we repeat this process to eventually get $_5K_{30}$ (20 years later, as this is coming from the given $_5K_{10}$)

and since we know the population has been growing at a fixed rate r for a long time, we know its a stable population and each age group population grows at the same rate as the population.

so the value we got for $_5K_{30}$ after 20 years can help us solve for the "current" value.

$_5K_{30}^{20 years later} = {_5}K_{30} \times e^{r\times 20}$ 

so what we want comes from calculation $\frac{_5K_{30}^{20 years later}}{e^{r\times20}}$ 