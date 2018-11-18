Renee Senining

# Extra Credit

1. First, I thought of the # of cars still on the road at each year as the "survivorship" with the radix being size 1,500,000 for both A and B.

   1. First, I thought of the # of cars still on the road at each year as the "survivorship" with the radix being size 1,500,000 for both A and B.

   2. Then, I calculated the number of "car deaths" $_nd_x$ by simply subtracting the number of cars in year that have been on the road for x+n  years from the number of cars that have been on the road for x years.

   3. Then, I calculated the probability the car breaks down with n years after being around for x years $_nq_x =  \frac{_nd_x}{l_x}$ and $_ np_x = 1 - {_n}q_x$ 

   4. I calculated the "car years lived/on the road" the way I calculate $_nL_x = \frac{n}{2}\times(l_x + l_{x+n})$ except for the last interval...

   5. For year 2016 $l_x$ , I made it equal to $0.3 \times 1500000$ for both A and B, and using the fact that the average time remaining after 16 years on the road for each car was 5, I found ${_\infin}L_{16} = 0 + l_{16}\times5$  

   6. Then I calculated $T_x$ and $e_x$ as usual and got the following tables

      | x    | $l_x^A$   | $_nd_x^A$ | $_nq_x^A$ | $_np_x^A$ | $_nL_x^A$ | $T_x^A$    | $e_x^A$      |
      | ---- | --------- | --------- | --------- | --------- | --------- | ---------- | ------------ |
      | 0    | 1,500,000 | 203,066   | 0.13538   | 0.86562   | 5,593,868 | 17,697,229 | **11.79815** |
      | 4    | 1,296,934 | 298,291   | 0.23000   | 0.77000   | 4,591,154 | 12,103,361 | 9.33229      |
      | 8    | 998,643   | 254,318   | 0.25466   | 0.74534   | 1,742,968 | 7,512,207  | 7.52241      |
      | 10   | 744,325   | 95,116    | 0.12779   | 0.87221   | 1,393,534 | 5,769,239  | 7.75097      |
      | 12   | 649,209   | 135,961   | 0.20943   | 0.79057   | 1,162,457 | 4,375,705  | 6.74006      |
      | 14   | 513,248   | 63,248    | 0.12323   | 0.87677   | 963,248   | 3,213,248  | 6.26061      |
      | 16   | 450,000   | 450,000   | 1         | 0         | 2,250,000 | 2,250,000  | 5.00000      |

      | x    | $l_x^B$   | $_nd_x^B$ | $_nq_x^B$ | $_np_x^B$ | $_nL_x^B$ | $T_x^B$    | $e_x^B$      |
      | ---- | --------- | --------- | --------- | --------- | --------- | ---------- | ------------ |
      | 0    | 1,500,000 | 178,773   | 0.11918   | 0.88082   | 5,642,454 | 17,396,747 | **11.59783** |
      | 4    | 1,321,227 | 375,982   | 0.28457   | 0.71543   | 4,532,944 | 11,754,293 | 8.89650      |
      | 8    | 945,245   | 258,691   | 0.27368   | 0.72632   | 1,631,799 | 7,221,349  | 7.63966      |
      | 10   | 686,554   | 84,601    | 0.12323   | 0.87677   | 1,288,507 | 5,589,550  | 8.14146      |
      | 12   | 601,953   | 102,408   | 0.17013   | 0.82987   | 1,101,498 | 4,301,043  | 7.14515      |
      | 14   | 499,545   | 49,545    | 0.09918   | 0.90082   | 949,545   | 4,199,545  | 6.40492      |
      | 16   | 450,000   | 450,000   | 1         | 0         | 2,250,000 | 2,250,000  | 5.00000      |

   Seems like Company A's cars are more reliable on average number of years a new car will last before it breaks down.

2. First, I need $_3p_6$, which I can get from obtaining $_2p_6$ from getting the square root of $_4p_4$  and $_1p_8$ from the square root of $_2p_8$ (Because we can assume $_1p_x = _1p_{x+1} = . . . {_1}p_{x+n-1}= (_np_x)^{\frac{1}{n}}$)

   Then, I use $_3q_6 = 1 - {_3}p_6 = 1 - {_2}p_6{_1}p_8 $

   From 1, I got that A has more reliable cars, so using A's probabilities...

   $1 - (0.77)^{\frac{1}{2}}\times(0.74534)^{\frac{1}{2}}$ = **0.24243** 

   For Company B: **0.27914 ** 

   

3. I made the assumption that $e_{16} = 5$ and $l_{16} = {_\infin}d_{16}$  

   | x    | $l_x^C$   | $_nd_x^C$ | $_nq_x^C$ | $_np_x^C$ | $_nL_x^C$ | $T_x^C$    | $e_x^C$      |
   | ---- | --------- | --------- | --------- | --------- | --------- | ---------- | ------------ |
   | 0    | 1,500,000 | 37,668    | 0.02511   | 0.97489   | 5,924,664 | 19,977,467 | **13.31831** |
   | 4    | 1,462,332 | 251,239   | 0.17181   | 0.82819   | 5,346,850 | 14,052,803 | 9.60986      |
   | 8    | 1,211,093 | 332,871   | 0.27485   | 0.72515   | 2,089,315 | 8,705,953  | 7.18851      |
   | 10   | 878,222   | 153,267   | 0.17452   | 0.82548   | 1,603,177 | 6,616,638  | 7.53413      |
   | 12   | 724,955   | 125,630   | 0.17329   | 0.82671   | 1,324,280 | 5,013,461  | 6.91555      |
   | 14   | 599,325   | 84,349    | 0.14704   | 0.85926   | 1,114,301 | 3,689,181  | 6.15556      |
   | 16   | 514,976   | 514,976   | 1         | 0         | 2,574,880 | 2,574,880  | 5.0000       |

   

Company C would do well, better than both A and B.