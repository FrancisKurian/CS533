(* ::Package:: *)

(*                           DFcriticalValues.m                            *)

  NormalRV:= Module[{x, y},
                     y = Random[];
                     x = Quantile[NormalDistribution[0, 1], y];
                    {x, y}]

  NormalRVgraph:= 
     Module[{x, y}, 
            y = Random[];
            Off[NSolve::ifun];
            x = NSolve[CDF[NormalDistribution[0, 1], a] == y, a][[1, 1, 2]];
            Show[Plot[CDF[NormalDistribution[0, 1], w], {w, -4, 4}, 
                      PlotRange -> {{-4, 4}, {0, 1}}], 
                 Graphics[{RGBColor[1, 0, 0, 1], 
                           Dashing[{0.005, 0.005, 0.005}], 
                           Line[{{0, y}, {x, y}, {x, 0}}]}], 
                 Graphics[Text[StyleForm["x = ", FontSlant -> Plain, 
                               FontFamily -> "TimesNewRoman", 16], 
                               {2, 0.5}]],
                 Graphics[Text[StyleForm[x, FontSlant -> Plain, 
                               FontFamily -> "TimesNewRoman", 16], 
                               {2.2, 0.5}, {-1, 0}]]]]

(*  The function NormalRVs[n] generates a sequence of 'n' indepentdent     *)
(*  normal random variables with mean zero and variance 2.  The variance   *)
(*  can be changed directly in the function definition.                    *)

    NormalRVs[n_, \[Sigma]_]:= Table[Random[NormalDistribution[0, \[Sigma]]], 
                                     {i, 1, n}]

(*  The function PriceSequence[] generates an autoregressive sequence of   *)
(*  length 'n' with adjustment rate 'c', initial price p0, equilibrium     *)
(*  price pStar.                                                           *)

    PriceSequence[c_, p0_, pStar_, n_, \[Sigma]_]:= 
       Module[{rv, p},
              rv = NormalRVs[n, \[Sigma]];
              p = {p0};
              Label[start];
              p = Join[p, {p[[-1]] + c * (pStar - p[[-1]]) + 
                           rv[[Length[p]]]}];
              If[Length[p] < n, Goto[start]];
              Return[p]]

    PriceSequenceGraph[c_, p0_, pStar_, n_, \[Sigma]_]:= 
        Module[{s, r, o}, 
                s = PriceSequence[c, p0, pStar, n, \[Sigma]]; 
                r = If[Min[s] > 50, 
                       {{0, n}, {10 * Floor[0.1 * Min[s]], 
                                  10 * Ceiling[0.1 * Max[s]]}},
                       {{0, n}, {0.05 * Floor[20 * Min[s]], 
                                  0.05 * Ceiling[20 * Max[s]]}}];
                Show[ListPlot[s, Joined -> True, PlotRange -> r, 
                              AxesOrigin -> {0, 0.05 * Floor[20 * Min[s]]}], 
                     Graphics[{Dashing[{0.01, 0.01, 0.01}], 
                               Line[{{0, pStar}, {n, pStar}}]}]]]













(*  Simulate a linear sequence with i.i.d. errors and then estimate        *)
(*  model parameters.                                                      *)

    f[x_, a_, b_, \[Sigma]_]:= a + b * x + 
                               Random[NormalDistribution[0, \[Sigma]]]

    bEsttStatLinearTest[a_, b_, n_, \[Sigma]_]:= 
       Module[{d, r, est, bHat, se, t, pHat},
              d = Table[f[i, a, b, \[Sigma]], {i, 1, n}];
              r = LinearModelFit[d, {x}, x];
            est = r["ParameterTable"];
           bHat = r["ParameterTable"][[1, 1, 3, 2]];
             se = r["ParameterTable"][[1, 1, 3, 3]];
              t = bHat/se; 
           Return[{d, r, est, bHat, se, t}]]

    bEsttStatLinear[a_, b_, n_, \[Sigma]_]:= 
       Module[{d, r, est, bHat, se, t, p},
              d = Table[f[i, a, b, \[Sigma]], {i, 1, n}];
              r = LinearModelFit[d, {x}, x];
            est = r["ParameterTable"];
           bHat = r["ParameterTable"][[1, 1, 3, 2]];
             se = r["ParameterTable"][[1, 1, 3, 3]];
              p = r["ParameterTable"][[1, 1, 3, 5]];
              t = bHat/se; 
           Return[{bHat, t, p}]]

    bEsttStatLinearGraph[a_, b_, n_, \[Sigma]_]:= 
       Module[{d, r, est, aHat, bHat, se, t, pHat, ra, o, g},
              d = Table[f[i, a, b, \[Sigma]], {i, 1, n}];
              r = LinearModelFit[d, {x}, x];
            est = r["ParameterTable"];
           aHat = r["ParameterTable"][[1, 1, 2, 2]];
           bHat = r["ParameterTable"][[1, 1, 3, 2]];
             se = r["ParameterTable"][[1, 1, 3, 3]];
              t = bHat/se; 
             ra = If[Min[d] < 50, 
                     {{0, n * 1.005}, {2 * Floor[0.5 * Min[d]], 
                                2 * Ceiling[0.5 * Max[d]]}},
                     {{0, n * 1.005}, {10 * Floor[0.1 * Min[d]], 
                                10 * Ceiling[0.1 * Max[d]]}}];
              o = If[Min[d] < 50, 2 * Floor[0.5 * Min[d]],
                                   10 * Floor[0.1 * Min[d]]];
              g = Show[ListPlot[d, Joined -> False, PlotRange -> ra, 
                                AxesLabel -> {x, y},
                                AxesOrigin -> {0, o},
                                ImageSize -> Medium,
                                AspectRatio -> 0.4],
                        Graphics[{Dashing[{0.01, 0.01, 0.01}], 
                                  Line[{{0, aHat}, {n, aHat + bHat * n}}]}]];
              Return[{est, g}]]

(*  The function cEsttStat reports the estimate of 'c'.                    *)

    cEsttStat[c_, p0_, pStar_, n_, \[Sigma]_, const_]:= 
       Module[{pr, d, r, est, a, se, t, pHat},
            pr = PriceSequence[c, p0, pStar, n, \[Sigma]];
            (* Dickey-Fuller estimates *) 
            d = Transpose[{-Drop[pr, -1], Drop[pr, 1] - Drop[pr, -1]}];
            r = If[const == 1, LinearModelFit[d, {x}, x], 
                               LinearModelFit[d, {x}, x, 
                                              IncludeConstantBasis -> False]];
            est = r["ParameterTable"];
            a = If[const == 1, r["ParameterTable"][[1, 1, 3, 2]], 
                               r["ParameterTable"][[1, 1, 2, 2]]];
            se = If[const == 1, r["ParameterTable"][[1, 1, 3, 3]], 
                                r["ParameterTable"][[1, 1, 2, 3]]];
            pHat = If[const == 1, est[[1, 1, 2, 2]]/(1 - est[[1, 1, 3, 2]])];
            t = a/se; 
            If[const == 1, Return[{a, t, pHat}], Return[{a, t}]]]

(*  The next two function provide feedback on intermediate stages of the   *)
(*  calculations in cEsttStat.                                             *)

    cEsttStatTest1[c_, p0_, pStar_, n_, \[Sigma]_, const_]:= 
       Module[{pr, d, r, est, a, se, t, pHat},
            pr = PriceSequence[c, p0, pStar, n, \[Sigma]];
            (* Dickey-Fuller estimates *) 
            d = Transpose[{-Drop[pr, -1], Drop[pr, 1] - Drop[pr, -1]}];
            r = If[const == 1, LinearModelFit[d, {x}, x], 
                               LinearModelFit[d, {x}, x, 
                                              IncludeConstantBasis -> False]];
            est = r["ParameterTable"];
            est]

    cEsttStatTest2[c_, p0_, pStar_, n_, \[Sigma]_, const_]:= 
       Module[{pr, d, r, est, a, se, t, pHat},
            pr = PriceSequence[c, p0, pStar, n, \[Sigma]];
            (* Dickey-Fuller estimates *) 
            d = Transpose[{-Drop[pr, -1], Drop[pr, 1] - Drop[pr, -1]}];
            r = If[const == 1, LinearModelFit[d, {x}, x], 
                               LinearModelFit[d, {x}, x, 
                                              IncludeConstantBasis -> False]];
            est = r["ParameterTable"];
            a = If[const == 1, r["ParameterTable"][[1, 1, 3, 2]], 
	                       r["ParameterTable"][[1, 1, 2, 2]]];
	    se = If[const == 1, r["ParameterTable"][[1, 1, 3, 3]], 
	                        r["ParameterTable"][[1, 1, 2, 3]]];
	    pHat = If[const == 1, est[[1, 1, 2, 2]]/est[[1, 1, 3, 2]]];
            t = a/se; 
            If[const == 1, Return[{est, a, t, est[[1, 1, 3, 4]], pHat}], 
                           Return[{a, t}]]]

(*  The function cEsttStatGraph returns the estimate of the adjustment     *)
(*  rate, the t-statistic for the adjustment rate, the value of pHat (if   *)
(*  const = 1), and a graph of the price sequence.                         *)

    cEsttStatGraph[c_, p0_, pStar_, n_, \[Sigma]_, const_]:= 
       Module[{pr, lhs, rhs, d, r, est, a, se, t, pHat, ra, g},
            Off[General::munfl];
            pr = PriceSequence[c, p0, pStar, n, \[Sigma]];
            (* Dickey-Fuller estimates *) 
            d = Transpose[{Drop[pr, -1], Drop[pr, 1]}];
            r = If[const == 1, LinearModelFit[d, {x}, x], 
                               LinearModelFit[d, {x}, x, 
                                              IncludeConstantBasis -> False]];
            est = r["ParameterTable"];
            a = If[const == 1, r["ParameterTable"][[1, 1, 3, 2]], 
                               r["ParameterTable"][[1, 1, 2, 2]]];
            se = If[const == 1, r["ParameterTable"][[1, 1, 3, 3]], 
                                r["ParameterTable"][[1, 1, 2, 3]]];
            pHat = If[const == 1, est[[1, 1, 2, 2]]/(1 - est[[1, 1, 3, 2]])];
            t = (1 - a)/se; 
            ra = If[Min[pr] > 50, 
                    {{0, n}, {10 * Floor[0.1 * Min[pr]], 
                              10 * Ceiling[0.1 * Max[pr]]}},
                    {{0, n}, {10 * Floor[0.1 * Min[pr]], 
                              10 * Ceiling[0.1 * Max[pr]]}}];
            g = Show[ListPlot[pr, Joined -> True, PlotRange -> ra, 
                              AxesOrigin -> {0, 10 * Floor[0.1 * Min[pr]]},
                              ImageSize -> Medium], 
                     Graphics[{Dashing[{0.01, 0.01, 0.01}], 
                               Line[{{0, pHat}, {n, pHat}}]}]];
            If[const == 1, Return[{1 - a, t, pHat, g}], Return[{1 - a, t, g}]]]



    cEsttStat[data_, const_]:= 
       Module[{d, r, est, a, se, t, pHat},
            (* Dickey-Fuller estimates *) 
            d = Transpose[{Drop[data, -1], Drop[data, 1]}];
	    r = If[const == 1, LinearModelFit[d, {x}, x], 
	                       LinearModelFit[d, {x}, x, 
	                                      IncludeConstantBasis -> False]];
	    est = r["ParameterTable"];
	    a = If[const == 1, r["ParameterTable"][[1, 1, 3, 2]], 
	                       r["ParameterTable"][[1, 1, 2, 2]]];
	    se = If[const == 1, r["ParameterTable"][[1, 1, 3, 3]], 
	                        r["ParameterTable"][[1, 1, 2, 3]]];
	    pHat = If[const == 1, est[[1, 1, 2, 2]]/(1 - est[[1, 1, 3, 2]])];
            t = (1 - a)/se; 
            If[const == 1, Return[{1 - a, t, pHat}], Return[{1 - a, t}]]]









(*  The function PriceSequenceE[] generates an autoregressive sequence of   *)
(*  length 'n' according to the same model as PriceSequenceE[] but the      *) 
(*  values of the prices are exponentiated.                                 *)

    PriceSequenceE[c_, p0_, pStar_, n_, \[Sigma]_]:= 
       Module[{rv, p},
              rv = NormalRVs[n, \[Sigma]];
              p = {p0};
              Label[start];
              p = Join[p, {p[[-1]] + c * (pStar - p[[-1]]) + 
                           rv[[Length[p]]]}];
              If[Length[p] < n, Goto[start]];
              Return[Map[Exp, p]]]

    PriceSequenceEGraph[c_, p0_, pStar_, n_, \[Sigma]_]:= 
        Module[{s, r, o}, 
                s = PriceSequenceE[c, p0, pStar, n, \[Sigma]]; 
                r = If[Min[s] > 50, 
                       {{0, n}, {10 * Floor[0.1 * Min[s]], 
                                  10 * Ceiling[0.1 * Max[s]]}},
                       {{0, n}, {0.05 * Floor[20 * Min[s]], 
                                  0.05 * Ceiling[20 * Max[s]]}}];
                Show[ListPlot[s, Joined -> True, PlotRange -> r, 
                              AxesOrigin -> {0, 0.05 * Floor[20 * Min[s]]}], 
                     Graphics[{Dashing[{0.01, 0.01, 0.01}], 
                               Line[{{0, E^pStar}, {n, E^pStar}}]}]]]

(*  The function cEsttStat reports the estimate of 'c' *)

    cEsttStatEL[c_, p0_, pStar_, n_, \[Sigma]_, const_]:= 
       Module[{prE, pr, lhs, rhs, d, r, est, a, se, t},
            prE = PriceSequenceE[c, p0, pStar, n, \[Sigma]];
            pr = Map[Log, prE];
            (* Dickey-Fuller estimates *) 
            lhs = Drop[pr, 1] - Drop[pr, -1];
            rhs = -(Drop[pr, -1] - pStar);
            d = Transpose[{rhs, lhs}];
            r = If[const == 1,  LinearModelFit[d, {x}, x], 
                   LinearModelFit[d, {x}, x, IncludeConstantBasis -> False]];
            est =  r["ParameterTable"];
            a = If[const == 1, r["ParameterTable"][[1, 1, 3, 2]], 
                   r["ParameterTable"][[1, 1, 2, 2]]];
            se = If[const == 1, r["ParameterTable"][[1, 1, 3, 3]], 
                    r["ParameterTable"][[1, 1, 2, 3]]];
            t = a/se;
            Return[{1 - a, t}]]

(*  The function cEsttStatE reports the estimate of 'c' *)

    EcEsttStat[c_, p0_, pStar_, n_, \[Sigma]_, const_]:= 
       Module[{prE, pr, lhs, rhs, d, r, est, a, se, t},
            pr = PriceSequenceE[c, p0, pStar, n, \[Sigma]];
            (* Dickey-Fuller estimates *) 
            lhs = Drop[pr, 1] - Drop[pr, -1];
            rhs = -(Drop[pr, -1] - pStar);
            d = Transpose[{rhs, lhs}];
            r = If[const == 1,  LinearModelFit[d, {x}, x], 
                   LinearModelFit[d, {x}, x, IncludeConstantBasis -> False]];
            est =  r["ParameterTable"];
            a = If[const == 1, r["ParameterTable"][[1, 1, 3, 2]], 
                   r["ParameterTable"][[1, 1, 2, 2]]];
            se = If[const == 1, r["ParameterTable"][[1, 1, 3, 3]], 
                    r["ParameterTable"][[1, 1, 2, 3]]];
            t = a/se;
            Return[{1 - a, t}]]
