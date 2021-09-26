(*                               ARlags.m                                  *)
(*                  Price adjustment process simulation                    *)

(*   PriceSequencesLags is a mean reverting adjustment process that also   *)
(*   includes lagged price changes.                                        *) 

     NormalRVs[n_, σ_]:= Table[Random[NormalDistribution[0, σ]], 
                                      {i, 1, n}]

     PriceSequenceLags[b_, c_, p0_, pStar_, n_, σ_]:= 
       Module[{rv, p, start},
              rv = NormalRVs[n, σ];
              p = {p0};
              Label[start];
              p = Join[p, {p[[-1]] + c * (pStar - p[[-1]]) + 
                           Sum[If[Length[p] > i, 
                                  b[[i]] * (p[[-i]] - p[[-i-1]]), 0], 
                               {i, 1, 9}] + 
                           rv[[Length[p]]]}];
              If[Length[p] < n, Goto[start]];
              p = Table[{i, p[[i + 1]]}, {i, 0, n - 1}];
              Return[p]]

    PriceSequenceGraphLags[b_, c_, p0_, pStar_, n_, σ_, join_]:= 
        Module[{s, r1, r, o}, 
                s = PriceSequenceLags[b, c, p0, pStar, n, σ]; 
                r1 = Transpose[s][[2]];
                r = If[Min[r1] > 50, 
                       {{0, n + 1}, {0.1 * Floor[10 * Min[0.98 * r1]], 
                                     0.1 * Ceiling[10 * Max[1.02 * r1]]}},
                       If[Min[r1] > 0, 
                           {{0, n + 1}, {0.05 * Floor[20 * Min[0.98 * r1]], 
                                          0.05 * Ceiling[20 * Max[1.025 * r1]]}},
                           {{0, n + 1}, {0.05 * Floor[20 * Min[1.02 * r1]],
                                          0.05 * Ceiling[20 * Max[1.02 * r1]]}}]];
                Show[ListPlot[s, Joined -> join, PlotRange -> r, 
                              AxesOrigin -> {r[[1, 1]], r[[2, 1]]}, 
                              PlotStyle -> PointSize[0.01]], 
                     Graphics[{Dashing[{0.01, 0.01, 0.01}], 
                               Line[{{0, pStar}, {n, pStar}}]}]]]

(*                  Estimate the adjustment rate and lags                  *)

    bcn[data_]:= 
      Module[{p, lhs, rhs, 
              dpL1, dpL2, dpL3, dpL4, dpL5, dpL6, dpL7, dpL8, dpL9, c},
             p = data;
             lhs = Drop[p, 1] - Drop[p, -1];
             rhs = - Drop[p, -1]; 
             dpL1 = Join[{0}, Drop[lhs, -1]];
             dpL2 = Join[{0, 0}, Drop[lhs, -2]];
             dpL3 = Join[{0, 0, 0}, Drop[lhs, -3]];
             dpL4 = Join[{0, 0, 0, 0}, Drop[lhs, -4]];
             dpL5 = Join[{0, 0, 0, 0, 0}, Drop[lhs, -5]];
             dpL6 = Join[{0, 0, 0, 0, 0, 0}, Drop[lhs, -6]];
             dpL7 = Join[{0, 0, 0, 0, 0, 0, 0}, Drop[lhs, -7]];
             dpL8 = Join[{0, 0, 0, 0, 0, 0, 0, 0}, Drop[lhs, -8]];
             dpL9 = Join[{0, 0, 0, 0, 0, 0, 0, 0, 0}, Drop[lhs, -9]];
             {rhs, dpL1, dpL2, dpL3, dpL4, dpL5, 
                   dpL6, dpL7, dpL8, dpL9, lhs}]

(*   The next function is similar to the function bcn[] but it drops the   *)
(*   terms when the lags would go back before the start of the sequence.   *)
(*   In bcn[] those lags are arbitrarily assigned a value of zero.        *)

     bcn2[data_, L_]:= 
       Module[{p, lhs, rhs, 
               dpL1, dpL2, dpL3, dpL4, dpL5, dpL6, dpL7, dpL8, dpL9, c},
              p = data;
              lhs = Drop[p, 1] - Drop[p, -1];
              rhs = - Drop[Drop[p, -1], L]; 
              dpL1 = Drop[Drop[lhs, -1], L - 1];
              dpL2 = Drop[Drop[lhs, -2], L - 2];
              dpL3 = Drop[Drop[lhs, -3], L - 3];
              dpL4 = Drop[Drop[lhs, -4], L - 4];
              dpL5 = Drop[Drop[lhs, -5], L - 5];
              dpL6 = Drop[Drop[lhs, -6], L - 6];
              dpL7 = Drop[Drop[lhs, -7], L - 7];
              dpL8 = Drop[Drop[lhs, -8], Max[L - 8, 0]];
              dpL9 = Drop[Drop[lhs, -9], Max[L - 9, 0]];
              lhs = Drop[lhs, L];
              {rhs, dpL1, dpL2, dpL3, dpL4, dpL5, 
                    dpL6, dpL7, dpL8, dpL9, lhs}]

(*   Test for the presence of lagged adjustment.  (Also compare pages      *)
(*   226-227 in Applied Econometric Time Series by Walter Enders.)         *)

     LagTestn[data_, n_] := 
       Module[{L, v1, v2, v, w1, w2, w, r, m1, b, c, pHat, s1, m2, s2}, 
              v1 = Take[bcn[data], n + 1];
              v2 = bcn[data][[-1]];
               v = Transpose[Join[v1, {v2}]];
              w1 = Take[bcn2[data, n], n + 1];
              w2 = bcn2[data, n][[-1]];
               w = Transpose[Join[w1, {w2}]];
               r = Take[{a, b1, b2, b3, b4, b5, b6, b7, b8, b9}, n + 1];
              m1 = LinearModelFit[v, r, r];
               b = m1["ParameterTable"][[1, 1, 2, 2]];
               c = m1["ParameterTable"][[1, 1, 3, 2]];
              pHat = b/c;
              s1 = Flatten[{m1[{"ParameterTable", "RSquared"}], pHat}];
              m2 = LinearModelFit[w, r, r];
              s2 = Flatten[{m2[{"ParameterTable", "RSquared"}], pHat}];
              Return[s2]]

(*                  Estimate the adjustment rate and lags                  *)

    bc[data_]:= 
      Module[{p, lhs, rhs, 
              dpL1, dpL2, dpL3, dpL4, dpL5, dpL6, dpL7, dpL8, dpL9, c},
             p = data;
             lhs = Drop[p, 1] - Drop[p, -1];
             rhs = - Drop[p, -1]; 
             dpL1 = Join[{0}, Drop[lhs, -1]];
             dpL2 = Join[{0, 0}, Drop[lhs, -2]];
             dpL3 = Join[{0, 0, 0}, Drop[lhs, -3]];
             dpL4 = Join[{0, 0, 0, 0}, Drop[lhs, -4]];
             dpL5 = Join[{0, 0, 0, 0, 0}, Drop[lhs, -5]];
             dpL6 = Join[{0, 0, 0, 0, 0, 0}, Drop[lhs, -6]];
             dpL7 = Join[{0, 0, 0, 0, 0, 0, 0}, Drop[lhs, -7]];
             dpL8 = Join[{0, 0, 0, 0, 0, 0, 0, 0}, Drop[lhs, -8]];
             dpL9 = Join[{0, 0, 0, 0, 0, 0, 0, 0, 0}, Drop[lhs, -9]];
             {rhs, dpL1, dpL2, dpL3, dpL4, dpL5, dpL6, dpL7, dpL8, dpL9, lhs}]

(*   Test for the presence of lagged adjustment.  (Also compare pages      *)
(*   226-227 in Applied Econometric Time Series by Walter Enders.)         *)

     LagTest[data_, L_]:= 
       Module[{v0, v1, v2, v, r, m, s, b, c, pHat}, 
              v0 = bc[data][[1]];
              v1 = Pick[Drop[Drop[bc[data], 1], -1], L, 1];
              v2 = bc[data][[-1]];
               v = Transpose[Join[{v0}, v1, {v2}]];
               r = Pick[{b1, b2, b3, b4, b5, b6, b7, b8, b9}, L, 1];
               r = Flatten[Join[{a}, r]];
               m = LinearModelFit[v, r, r];
               b = m["ParameterTable"][[1, 1, 2, 2]];
               c = m["ParameterTable"][[1, 1, 3, 2]];
              pHat = b/c;
               s = Flatten[{m[{"ParameterTable", "RSquared"}], pHat}];
               Return[s]]

(*   The next function collects parameter estimates from large numbers    *)
(*   of simulations.                                                      *)
(*   It returns estimates of the intercept 'b', the adjustment rate       *)
(*   'c', the lags 'l', the t-statistics for each of these estimates,     *)
(*   and the implied estimate of the level pHat that the series           *)
(*   fluctuates around.                                                   *)

     LagTestnDistribution[sim_, L_] := 
       Module[{N, s, est, b, c, l, tb, tc, tl, pHat, pl}, 
              N = Length[sim];
              s = Table[Transpose[sim[[i]]][[2]], {i, 1, N}]; 
              est = Table[LagTestn[s[[i]], L][[1, 1]], {i, 1, N}];
              b = Table[est[[i, 1, 2, 2]], {i, 1, N}];
              c = Table[est[[i, 1, 3, 2]], {i, 1, N}];
              l = Table[Table[est[[i, 1, 3 + j, 2]], {i, 1, N}], {j, 1, L}]; 
              tb = Table[est[[i, 1, 2, 4]], {i, 1, N}];
              tc = Table[est[[i, 1, 3, 4]], {i, 1, N}];
              tl = Table[Table[est[[i, 1, 3 + j, 4]], {i, 1, N}], {j, 1, L}];
              pl = Table[Table[est[[i, 1, 3 + j, 5]], {i, 1, N}], {j, 1, L}];
              pHat =  Table[b[[i]]/c[[i]], {i, 1, N}]; 
              Return[{b, c, l, tb, tc, tl, pHat}]]
