(* ::Package:: *)

(*                                Spurious.m                               *)

  NormalRV:= Module[{x, y},
                     y = Random[];
                     x = Quantile[NormalDistribution[0, 1], y];
                    {x, y}]

  NormalRVgraph:= 
     Module[{x, y}, 
            y = Random[];
            x = NSolve[CDF[NormalDistribution[0, 1], a] == y, a][[1, 1, 2]];
            Show[Plot[CDF[NormalDistribution[0, 1], w], {w, -4, 4}, 
                      PlotRange -> {{-4, 4}, {0, 1}}], 
                
                 Graphics[Text[StyleForm["x = ", FontSlant -> Plain, 
                               FontFamily -> "TimesNewRoman", 16], 
                               {2, 0.5}]],
                 Graphics[Text[StyleForm[x, FontSlant -> Plain, 
                               FontFamily -> "TimesNewRoman", 16], 
                               {2.2, 0.5}, {-1, 0}]]]]

(*  The function NormalRVs[n] generates a sequence of 'n' independent      *)
(*  normal random variables with mean zero and variance 2.  The variance   *)
(*  can be changed directly in the function definition.                    *)

    NormalRVs[n_]:= Table[Random[NormalDistribution[0, 2]], {i, 1, n}]

(*  The function PriceSequence[] generates an autoregressive sequence of   *)
(*  length 'n' with adjustment rate 'c', initial price p0, equilibrium     *)
(*  price pStar.                                                           *)

    PriceSequence[c_, p0_, pStar_, n_]:= 
       Module[{rv, p},
              rv = NormalRVs[n];
              p = {p0};
              Label[start];
              p = Join[p, {p[[-1]] + c * (pStar - p[[-1]]) + 
                           rv[[Length[p]]]}];
              If[Length[p] < n, Goto[start]];
              Return[p]]

(*  The function SpuriousEst reports the estimate of returns the two       *)
(*  sequences, the paired values from them ('d'), and the estimates from   *)
(*  the regression of one sequence on the other.                           *)

    SpuriousEst[c_, p0_, pStar_, n_, const_]:= 
       Module[{pr, lhs, rhs, d, r, est, int, slope, lhsRange, rhsRange, p},
            lhs = PriceSequence[c, p0, pStar, n];
            lhsRange = {5 * Floor[0.2 * Min[lhs]], 
                        5 * Ceiling[0.2 * Max[lhs]]};
            rhs = PriceSequence[c, p0, pStar, n];
            rhsRange = {5 * Floor[0.2 * Min[rhs]], 
                        5 * Ceiling[0.2 * Max[rhs]]};
            d = Transpose[{rhs, lhs}];
            r = If[const == 1,  LinearModelFit[d, {x}, x], 
                   LinearModelFit[d, {x}, x, IncludeConstantBasis -> False],
                                             RegressionReport -> {DurbinWatsonD}];
            est = r["ParameterTable"]; 
            p = r["ANOVATable"];
            Return[{lhs, rhs, d, est, r[[1, 2]], lhsRange, rhsRange, p[[1, 1, 2, 6]]}]]

(*  The function f[d] determines the endpoints of the regression line for  *) 
(*  the random walks that comprise the variables 'x' and 'y' in the        *)
(*  scatter plot.                                                          *)

    f[d_]:= Module[{a0, a1, T}, 
               a0 = LinearModelFit[d, {y}, y][[1, 2, 1]];
               a1 = LinearModelFit[d, {y}, y][[1, 2, 2]];
                T = Length[d];
               {{0, a0}, {T, a0 + T * a1}}]

(*  What is the Durbin-Watson statistic for each random walk separately?   *)

    Durbin[data_]:= Module[{a0, a1, T, e, d, d0}, 
		    a0 = LinearModelFit[data, {y}, y][[1, 2, 1]]; 
		    a1 = LinearModelFit[data, {y}, y][[1, 2, 2]];
             T = Length[data];
             e = Table[a0 + a1 * i - data[[i]], {i, 1, T}];
             d = Sum[(e[[t]] - e[[t - 1]])^2, {t, 2, T}]/Sum[e[[t]]^2, {t, 1, T}];
             d]

(*  This function determines the percentage of regressions that produce    *)
(*  a significant relationship between x and y.                            *)

    SpuriousSignificance[c_, p0_, pStar_, n_, const_, T_]:=
       Module[{t, s, u, L},
               t = ParallelTable[Quiet[SpuriousEst[c, p0, pStar, n, 1]][[8]], {T}];
               s = Sort[t] - 0.05;
               u = Sign[s];
               L = Length[Position[u, -1]];
               N[L/T]]

(*  Find the relation between the law of large numbers and the tendency    *)
(*  of a random walk to drift around.                                      *)
