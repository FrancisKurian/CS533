    (*                          Market Excess Demand                       *)

    (*                           Equilibrium Prices                        *)

    (*   This function calculates the equilibrium prices by subdivision    *)
    (*   of the price set.  The prices are first renormalized to use       *)
    (*   the unit simplex (px > 0, py > 0, px + py = 1) and then, once     *)
    (*   the equilibrium price is determined, the price is normalized      *)
    (*   again so that px = 1.  N = 256 should be an adequate number of    *)
    (*   subdivisions of the price set to get an accurate value for the    *)
    (*   equilibrium price.                                                *)

         (*  'sign' -  Determine whether ZY is positive or negative at     *)
         (*            each value on the subdivision.                      *)
         (*  'diff' -  Form an array of differences to determine cross-    *)
         (*            over points where ZY changes from negative to       *)
         (*            positive (or vice versa).  Each value in this       *)
         (*            table will be either 0 or 2.  Crossovers occur      *)
         (*            between the positions where a 2 appears and the     *)
         (*            next grid point to the right of that value.         *)
         (*  'l'    -  Find the left endpoint of a crossover interval.     *)
         (*  'lft'  -  Convert 'l' to a price vector of the form           *)
         (*            (p, 1 - p).                                         *)
         (*  'r'    -  Find the right endpoint of a crossover interval.    *)
         (*  'rht'  -  Convert 'r' to a price vector of the form           *)
         (*            (p, 1 - p).                                         *)
         (*  'fac'  -  Find the factor that we be used to weight the       *)
         (*            left (l) and right (r) endpoints in the linear      *)
         (*            interpolation that is used to determine the         *)
         (*            equilibrium price.                                  *)
         (*  'eqpr' -  Find the equilibrium price on the unit simplex.     *)
         (*            This is converted to the normalization used in      *)
         (*            the experiment by the calculation                   *)
         (*            p = (1 - eqpr)/eqpr.                                *)

         EqPrice[alpha_, gamma_, beta_, sB_, N_]:= 
              Module[{sign, diff, l, lft, r, rht, fac}, 
                     sign = Table[Sign[ZY[alpha, gamma, beta, sB, 
                                   {0.0001 + (0.9998 (i - 1))/(N - 1), 
                                    0.9999 - (0.9998 (i - 1))/(N - 1)}]],
                                  {i, N}];
                     diff = Abs[Drop[sign, 1] - Drop[sign, -1]];
                     l = 0.0001 + 0.9998 
                         * (Flatten[Position[diff, 2]][[1]] - 1)/(N - 1);
                     lft = {l, 1 - l};
                     r = l + 0.9998/(N - 1); 
                     rht = {r, 1 - r};
                     fac = Abs[ZY[alpha, gamma, beta, sB, lft]/
                              (ZY[alpha, gamma, beta, sB, lft] -
                               ZY[alpha, gamma, beta, sB, rht])];
                     eqpr = (1 - fac) * lft[[1]] + fac * rht[[1]];
                     Return[(1 - eqpr)/eqpr]]