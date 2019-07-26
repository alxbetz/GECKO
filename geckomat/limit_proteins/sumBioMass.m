%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,P,C,R,D,L] = sumBioMass(model)
% Calculates breakdown of biomass for the yeast model:
% X -> Biomass fraction without lipids [g/gDW]
% P -> Protein fraction [g/gDW]
% C -> Carbohydrate fraction [g/gDW]
% R -> RNA fraction [g/gDW]
% D -> DNA fraction [g/gDW]
% L -> Lipid fraction [g/gDW]
% Function adapted from SLIMEr: https://github.com/SysBioChalmers/SLIMEr
%
% Benjamin Sanchez. Last update: 2018-10-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,P,C,R,D,L] = sumBioMass(model)

%Components of biomass:
%       id         MW [g/mol]  class     name
% comps = {'s_0404'	89.09       'P'     % A     Alanine         ala
%          's_0542'	121.16      'P'     % C     Cysteine        cys
%          's_0432'   133.11      'P'     % D     Aspartic acid   asp
%          's_0748'   147.13      'P'     % E     Glutamic acid   glu
%          's_1314'   165.19      'P'     % F     Phenylalanine   phe
%          's_0757'   75.07       'P'     % G     Glycine         gly
%          's_0832'   155.15      'P'     % H     Histidine       his
%          's_0847'   131.17      'P'     % I     Isoleucine      ile
%          's_1099'   146.19      'P'     % K     Lysine          lys
%          's_1077'   131.17      'P'     % L     Leucine         leu
%          's_1148'   149.21      'P'     % M     Methionine      met
%          's_0430'   132.12      'P'     % N     Asparagine      asn
%          's_1379'   115.13      'P'     % P     Proline         pro
%          's_0747'   146.14      'P'     % Q     Glutamine       gln
%          's_0428'   174.2       'P'     % R     Arginine        arg
%          's_1428'   105.09      'P'     % S     Serine          ser
%          's_1491'   119.12      'P'     % T     Threonine       thr
%          's_1561'   117.15      'P'     % V     Valine          val
%          's_1527'   204.23      'P'     % W     Tryptophan      trp
%          's_1533'   181.19      'P'     % Y     Tyrosine        tyr
%          's_0001'	180.16      'C'     % (1->3)-beta-D-glucan
%          's_0004'	180.16      'C'     % (1->6)-beta-D-glucan
%          's_0509'   221.21      'C'     % chitin
%          's_0773'   180.16      'C'     % glycogen
%          's_1107'   180.16      'C'     % mannan
%          's_1520'   342.296 	'C'     % trehalose
%          's_0423'   347.22      'R'     % AMP
%          's_0526'   323.2       'R'     % CMP
%          's_0782'   363.22      'R'     % GMP
%          's_1545'   324.18      'R'     % UMP
%          's_0584'   331.22      'D'     % dAMP
%          's_0589'   307.2       'D'     % dCMP
%          's_0615'   345.21      'D'     % dGMP
%          's_0649'   322.21      'D'     % dTMP
%          's_3714'   852.83      'N'     % heme a
%          's_1405'   376.36      'N'     % riboflavin
%          's_1467'   96.06       'N'};   % sulphate





% 
% biomassMets = findMetsFromRxns(crem,{'Biomass_Chlamy_auto'});
% [MW, Ematrix]  = computeMW(crem,biomassMets);
% crem.metNames(find(ismember(crem.mets,biomassMets)))
%       id         MW [g/mol]  class     name
% comps = {'ac[c]'	59	'C'
% 'acaro[h]'	536	'N'
% 'adp[c]'	424	'N'
% 'alatrna[c]'	399	'P'
% 'anxan[u]'	584	'C'
% 'arab_L[c]'	150	'C'
% 'argtrna[c]'	619	'P'
% 'asntrna[c]'	442	'P'
% 'asptrna[c]'	442	'P'
% 'asqdca18111Z160[c]'	1077	'L'
% 'asqdca1819Z160[c]'	1077	'L'
% 'asqdca1829Z12Z160[c]'	1075	'L'
% 'asqdca1839Z12Z15Z160[c]'	1073	'L'
% 'asqdpa18111Z160[c]'	1079	'L'
% 'asqdpa1819Z160[c]'	1079	'L'
% 'asqdpa1829Z12Z160[c]'	1077	'L'
% 'asqdpa1839Z12Z15Z160[c]'	1075	'L'
% 'atp[c]'	503	'R'
% 'btn[c]'	243	'N'
% 'but[c]'	87	'L'
% 'caro[u]'	536	'N'
% 'chla[u]'	892	'N'
% 'chlb[u]'	906	'N'
% 'ctp[c]'	479	'N'
% 'cystrna[c]'	565	'P'
% 'datp[c]'	487	'D'
% 'dctp[c]'	463	'D'
% 'dgdg1819Z160[h]'	918	'L'
% 'dgdg1819Z1617Z[h]'	916	'L'
% 'dgdg1819Z1619Z[h]'	916	'L'
% 'dgdg1819Z1627Z10Z[h]'	914	'L'
% 'dgdg1819Z1634Z7Z10Z[h]'	912	'L'
% 'dgdg1819Z1637Z10Z13Z[h]'	912	'L'
% 'dgdg1829Z12Z160[h]'	916	'L'
% 'dgdg1829Z12Z1617Z[h]'	914	'L'
% 'dgdg1829Z12Z1619Z[h]'	914	'L'
% 'dgdg1829Z12Z1627Z10Z[h]'	912	'L'
% 'dgdg1829Z12Z1634Z7Z10Z[h]'	910	'L'
% 'dgdg1829Z12Z1637Z10Z13Z[h]'	910	'L'
% 'dgdg1839Z12Z15Z160[h]'	914	'L'
% 'dgdg1839Z12Z15Z1627Z10Z[h]'	910	'L'
% 'dgdg1839Z12Z15Z1634Z7Z10Z[h]'	908	'L'
% 'dgdg1839Z12Z15Z1637Z10Z13Z[h]'	908	'L'
% 'dgdg1839Z12Z15Z1644Z7Z10Z13Z[h]'	906	'L'
% 'dgtp[c]'	503	'D'
% 'dgts16018111Z[c]'	737	'L'
% 'dgts1601819Z[c]'	737	'L'
% 'dgts1601829Z12Z[c]'	735	'L'
% 'dgts1601835Z9Z12Z[c]'	733	'L'
% 'dgts1601845Z9Z12Z15Z[c]'	731	'L'
% 'dgts18111Z18111Z[c]'	763	'L'
% 'dgts18111Z1819Z[c]'	763	'L'
% 'dgts18111Z1829Z12Z[c]'	761	'L'
% 'dgts18111Z1835Z9Z12Z[c]'	759	'L'
% 'dgts18111Z1845Z9Z12Z15Z[c]'	757	'L'
% 'dgts1819Z18111Z[c]'	763	'L'
% 'dgts1819Z1819Z[c]'	763	'L'
% 'dgts1819Z1829Z12Z[c]'	761	'L'
% 'dgts1819Z1835Z9Z12Z[c]'	759	'L'
% 'dgts1819Z1845Z9Z12Z15Z[c]'	757	'L'
% 'dgts1829Z12Z18111Z[c]'	761	'L'
% 'dgts1829Z12Z1819Z[c]'	761	'L'
% 'dgts1829Z12Z1829Z12Z[c]'	759	'L'
% 'dgts1829Z12Z1835Z9Z12Z[c]'	757	'L'
% 'dgts1829Z12Z1845Z9Z12Z15Z[c]'	755	'L'
% 'dgts1839Z12Z15Z18111Z[c]'	759	'L'
% 'dgts1839Z12Z15Z1819Z[c]'	759	'L'
% 'dgts1839Z12Z15Z1835Z9Z12Z[c]'	755	'L'
% 'dgts1839Z12Z15Z1845Z9Z12Z15Z[c]'	753	'L'
% 'dttp[c]'	478	'D'
% 'fad[c]'	783	'N'
% 'gal[c]'	180	'C'
% 'glntrna[c]'	590	'P'
% 'glutrna[c]'	590	'P'
% 'glyc[c]'	92	'L'
% 'glytrna[c]'	385	'P'
% 'gthrd[c]'	305	'N'
% 'gtp[c]'	519	'R'
% 'h[c]'	1	'N'
% 'h2o[c]'	18	'N'
% 'histrna[c]'	465	'P'
% 'iletrna[c]'	575	'P'
% 'leutrna[c]'	575	'P'
% 'loroxan[u]'	584	'N'
% 'lut[u]'	568	'N'
% 'lystrna[c]'	457	'P'
% 'man[c]'	180	'C'
% 'mettrna[c]'	593	'P'
% 'mgdg1829Z12Z160[h]'	754	'L'
% 'mgdg1829Z12Z1617Z[h]'	752	'L'
% 'mgdg1829Z12Z1619Z[h]'	752	'L'
% 'mgdg1829Z12Z1627Z10Z[h]'	750	'L'
% 'mgdg1829Z12Z1634Z7Z10Z[h]'	748	'L'
% 'mgdg1829Z12Z1637Z10Z13Z[h]'	748	'L'
% 'mgdg1829Z12Z1644Z7Z10Z13Z[h]'	746	'L'
% 'mgdg1839Z12Z15Z160[h]'	752	'L'
% 'mgdg1839Z12Z15Z1627Z10Z[h]'	748	'L'
% 'mgdg1839Z12Z15Z1634Z7Z10Z[h]'	746	'L'
% 'mgdg1839Z12Z15Z1637Z10Z13Z[h]'	746	'L'
% 'mgdg1839Z12Z15Z1644Z7Z10Z13Z[h]'	744	'L'
% 'nad[c]'	662	'N'
% 'nadh[c]'	663	'N'
% 'nadp[c]'	740	'N'
% 'nadph[c]'	741	'N'
% 'neoxan[u]'	600	'N'
% 'pail18111Z160[c]'	836	'L'
% 'pail1819Z160[c]'	836	'L'
% 'pe1801835Z9Z12Z[c]'	741	'L'
% 'pe1801845Z9Z12Z15Z[c]'	739	'L'
% 'pe18111Z1835Z9Z12Z[c]'	739	'L'
% 'pe18111Z1845Z9Z12Z15Z[c]'	737	'L'
% 'pe1819Z1835Z9Z12Z[c]'	739	'L'
% 'pe1819Z1845Z9Z12Z15Z[c]'	737	'L'
% 'pe1829Z12Z1835Z9Z12Z[c]'	737	'L'
% 'pg18111Z160[h]'	748	'L'
% 'pg18111Z1613E[h]'	746	'L'
% 'pg1819Z160[h]'	748	'L'
% 'pg1819Z1613E[h]'	746	'L'
% 'pg1829Z12Z160[h]'	746	'L'
% 'pg1829Z12Z1613E[h]'	744	'L'
% 'pg1839Z12Z15Z160[h]'	744	'L'
% 'pg1839Z12Z15Z1613E[h]'	742	'L'
% 'phetrna[c]'	475	'P'
% 'pi[c]'	96	'N'
% 'ppa[c]'	73	'C'
% 'protrna[c]'	425	'P'
% 'rhodopsin[s]'	268	'N'
% 'sertrna[c]'	415	'P'
% 'sqdg160[h]'	793	'L'
% 'sqdg18111Z160[h]'	819	'L'
% 'sqdg1819Z160[h]'	819	'L'
% 'sqdg1829Z12Z160[h]'	817	'L'
% 'sqdg1839Z12Z15Z160[h]'	815	'L'
% 'starch300[h]'	48618	'C'
% 'tag16018111Z160[c]'	832	'L'
% 'tag16018111Z180[c]'	860	'L'
% 'tag16018111Z18111Z[c]'	858	'L'
% 'tag16018111Z1819Z[c]'	858	'L'
% 'tag16018111Z1835Z9Z12Z[c]'	854	'L'
% 'tag16018111Z1845Z9Z12Z15Z[c]'	852	'L'
% 'tag1601819Z160[c]'	832	'L'
% 'tag1601819Z180[c]'	860	'L'
% 'tag1601819Z18111Z[c]'	858	'L'
% 'tag1601819Z1819Z[c]'	858	'L'
% 'tag1601819Z1835Z9Z12Z[c]'	854	'L'
% 'tag1601819Z1845Z9Z12Z15Z[c]'	852	'L'
% 'tag1801819Z160[c]'	860	'L'
% 'tag1801819Z180[c]'	888	'L'
% 'tag1801819Z18111Z[c]'	886	'L'
% 'tag1801819Z1819Z[c]'	886	'L'
% 'tag1801819Z1835Z9Z12Z[c]'	882	'L'
% 'tag1801819Z1845Z9Z12Z15Z[c]'	880	'L'
% 'tag18111Z18111Z160[c]'	858	'L'
% 'tag18111Z18111Z180[c]'	886	'L'
% 'tag18111Z18111Z18111Z[c]'	884	'L'
% 'tag18111Z18111Z1819Z[c]'	884	'L'
% 'tag18111Z18111Z1835Z9Z12Z[c]'	880	'L'
% 'tag18111Z18111Z1845Z9Z12Z15Z[c]'	878	'L'
% 'tag18111Z1819Z160[c]'	858	'L'
% 'tag18111Z1819Z180[c]'	886	'L'
% 'tag18111Z1819Z18111Z[c]'	884	'L'
% 'tag18111Z1819Z1819Z[c]'	884	'L'
% 'tag18111Z1819Z1835Z9Z12Z[c]'	880	'L'
% 'tag18111Z1819Z1845Z9Z12Z15Z[c]'	878	'L'
% 'tag1819Z18111Z160[c]'	858	'L'
% 'tag1819Z18111Z180[c]'	886	'L'
% 'tag1819Z18111Z18111Z[c]'	884	'L'
% 'tag1819Z18111Z1819Z[c]'	884	'L'
% 'tag1819Z18111Z1835Z9Z12Z[c]'	880	'L'
% 'tag1819Z18111Z1845Z9Z12Z15Z[c]'	878	'L'
% 'tag1819Z1819Z160[c]'	858	'L'
% 'tag1819Z1819Z180[c]'	886	'L'
% 'tag1819Z1819Z18111Z[c]'	884	'L'
% 'tag1819Z1819Z1819Z[c]'	884	'L'
% 'tag1819Z1819Z1835Z9Z12Z[c]'	880	'L'
% 'tag1819Z1819Z1845Z9Z12Z15Z[c]'	878	'L'
% 'thmmp[c]'	343	'N'
% 'thrtrna[c]'	429	'P'
% 'trnaala[c]'	328	'P'
% 'trnaarg[c]'	462	'P'
% 'trnaasn[c]'	328	'P'
% 'trnaasp[c]'	328	'P'
% 'trnacys[c]'	462	'P'
% 'trnagln[c]'	462	'P'
% 'trnaglu[c]'	462	'P'
% 'trnagly[c]'	328	'P'
% 'trnahis[c]'	328	'P'
% 'trnaile[c]'	462	'P'
% 'trnaleu[c]'	462	'P'
% 'trnalys[c]'	328	'P'
% 'trnamet[c]'	462	'P'
% 'trnaphe[c]'	328	'P'
% 'trnapro[c]'	328	'P'
% 'trnaser[c]'	328	'P'
% 'trnathr[c]'	328	'P'
% 'trnatrp[c]'	462	'P'
% 'trnatyr[c]'	462	'P'
% 'trnaval[c]'	462	'P'
% 'trptrna[c]'	648	'P'
% 'tyrtrna[c]'	625	'P'
% 'utp[c]'	480	'N'
% 'valtrna[c]'	561	'P'
% 'vioxan[u]'	600	'N'
% 'zaxan[u]'	568	'N'};
comps = {'ac_c'	59	'C'
'acaro_h'	536	'N'
'adp_c'	424	'N'
'alatrna_c'	399	'P'
'anxan_u'	584	'C'
'arab_L_c'	150	'C'
'argtrna_c'	619	'P'
'asntrna_c'	442	'P'
'asptrna_c'	442	'P'
'asqdca18111Z160_c'	1077	'L'
'asqdca1819Z160_c'	1077	'L'
'asqdca1829Z12Z160_c'	1075	'L'
'asqdca1839Z12Z15Z160_c'	1073	'L'
'asqdpa18111Z160_c'	1079	'L'
'asqdpa1819Z160_c'	1079	'L'
'asqdpa1829Z12Z160_c'	1077	'L'
'asqdpa1839Z12Z15Z160_c'	1075	'L'
'atp_c'	503	'R'
'btn_c'	243	'N'
'but_c'	87	'L'
'caro_u'	536	'N'
'chla_u'	892	'N'
'chlb_u'	906	'N'
'ctp_c'	479	'N'
'cystrna_c'	565	'P'
'datp_c'	487	'D'
'dctp_c'	463	'D'
'dgdg1819Z160_h'	918	'L'
'dgdg1819Z1617Z_h'	916	'L'
'dgdg1819Z1619Z_h'	916	'L'
'dgdg1819Z1627Z10Z_h'	914	'L'
'dgdg1819Z1634Z7Z10Z_h'	912	'L'
'dgdg1819Z1637Z10Z13Z_h'	912	'L'
'dgdg1829Z12Z160_h'	916	'L'
'dgdg1829Z12Z1617Z_h'	914	'L'
'dgdg1829Z12Z1619Z_h'	914	'L'
'dgdg1829Z12Z1627Z10Z_h'	912	'L'
'dgdg1829Z12Z1634Z7Z10Z_h'	910	'L'
'dgdg1829Z12Z1637Z10Z13Z_h'	910	'L'
'dgdg1839Z12Z15Z160_h'	914	'L'
'dgdg1839Z12Z15Z1627Z10Z_h'	910	'L'
'dgdg1839Z12Z15Z1634Z7Z10Z_h'	908	'L'
'dgdg1839Z12Z15Z1637Z10Z13Z_h'	908	'L'
'dgdg1839Z12Z15Z1644Z7Z10Z13Z_h'	906	'L'
'dgtp_c'	503	'D'
'dgts16018111Z_c'	737	'L'
'dgts1601819Z_c'	737	'L'
'dgts1601829Z12Z_c'	735	'L'
'dgts1601835Z9Z12Z_c'	733	'L'
'dgts1601845Z9Z12Z15Z_c'	731	'L'
'dgts18111Z18111Z_c'	763	'L'
'dgts18111Z1819Z_c'	763	'L'
'dgts18111Z1829Z12Z_c'	761	'L'
'dgts18111Z1835Z9Z12Z_c'	759	'L'
'dgts18111Z1845Z9Z12Z15Z_c'	757	'L'
'dgts1819Z18111Z_c'	763	'L'
'dgts1819Z1819Z_c'	763	'L'
'dgts1819Z1829Z12Z_c'	761	'L'
'dgts1819Z1835Z9Z12Z_c'	759	'L'
'dgts1819Z1845Z9Z12Z15Z_c'	757	'L'
'dgts1829Z12Z18111Z_c'	761	'L'
'dgts1829Z12Z1819Z_c'	761	'L'
'dgts1829Z12Z1829Z12Z_c'	759	'L'
'dgts1829Z12Z1835Z9Z12Z_c'	757	'L'
'dgts1829Z12Z1845Z9Z12Z15Z_c'	755	'L'
'dgts1839Z12Z15Z18111Z_c'	759	'L'
'dgts1839Z12Z15Z1819Z_c'	759	'L'
'dgts1839Z12Z15Z1835Z9Z12Z_c'	755	'L'
'dgts1839Z12Z15Z1845Z9Z12Z15Z_c'	753	'L'
'dttp_c'	478	'D'
'fad_c'	783	'N'
'gal_c'	180	'C'
'glntrna_c'	590	'P'
'glutrna_c'	590	'P'
'glyc_c'	92	'L'
'glytrna_c'	385	'P'
'gthrd_c'	305	'N'
'gtp_c'	519	'R'
'h_c'	1	'N'
'h2o_c'	18	'N'
'histrna_c'	465	'P'
'iletrna_c'	575	'P'
'leutrna_c'	575	'P'
'loroxan_u'	584	'N'
'lut_u'	568	'N'
'lystrna_c'	457	'P'
'man_c'	180	'C'
'mettrna_c'	593	'P'
'mgdg1829Z12Z160_h'	754	'L'
'mgdg1829Z12Z1617Z_h'	752	'L'
'mgdg1829Z12Z1619Z_h'	752	'L'
'mgdg1829Z12Z1627Z10Z_h'	750	'L'
'mgdg1829Z12Z1634Z7Z10Z_h'	748	'L'
'mgdg1829Z12Z1637Z10Z13Z_h'	748	'L'
'mgdg1829Z12Z1644Z7Z10Z13Z_h'	746	'L'
'mgdg1839Z12Z15Z160_h'	752	'L'
'mgdg1839Z12Z15Z1627Z10Z_h'	748	'L'
'mgdg1839Z12Z15Z1634Z7Z10Z_h'	746	'L'
'mgdg1839Z12Z15Z1637Z10Z13Z_h'	746	'L'
'mgdg1839Z12Z15Z1644Z7Z10Z13Z_h'	744	'L'
'nad_c'	662	'N'
'nadh_c'	663	'N'
'nadp_c'	740	'N'
'nadph_c'	741	'N'
'neoxan_u'	600	'N'
'pail18111Z160_c'	836	'L'
'pail1819Z160_c'	836	'L'
'pe1801835Z9Z12Z_c'	741	'L'
'pe1801845Z9Z12Z15Z_c'	739	'L'
'pe18111Z1835Z9Z12Z_c'	739	'L'
'pe18111Z1845Z9Z12Z15Z_c'	737	'L'
'pe1819Z1835Z9Z12Z_c'	739	'L'
'pe1819Z1845Z9Z12Z15Z_c'	737	'L'
'pe1829Z12Z1835Z9Z12Z_c'	737	'L'
'pg18111Z160_h'	748	'L'
'pg18111Z1613E_h'	746	'L'
'pg1819Z160_h'	748	'L'
'pg1819Z1613E_h'	746	'L'
'pg1829Z12Z160_h'	746	'L'
'pg1829Z12Z1613E_h'	744	'L'
'pg1839Z12Z15Z160_h'	744	'L'
'pg1839Z12Z15Z1613E_h'	742	'L'
'phetrna_c'	475	'P'
'pi_c'	96	'N'
'ppa_c'	73	'C'
'protrna_c'	425	'P'
'rhodopsin_s'	268	'N'
'sertrna_c'	415	'P'
'sqdg160_h'	793	'L'
'sqdg18111Z160_h'	819	'L'
'sqdg1819Z160_h'	819	'L'
'sqdg1829Z12Z160_h'	817	'L'
'sqdg1839Z12Z15Z160_h'	815	'L'
'starch300_h'	48618	'C'
'tag16018111Z160_c'	832	'L'
'tag16018111Z180_c'	860	'L'
'tag16018111Z18111Z_c'	858	'L'
'tag16018111Z1819Z_c'	858	'L'
'tag16018111Z1835Z9Z12Z_c'	854	'L'
'tag16018111Z1845Z9Z12Z15Z_c'	852	'L'
'tag1601819Z160_c'	832	'L'
'tag1601819Z180_c'	860	'L'
'tag1601819Z18111Z_c'	858	'L'
'tag1601819Z1819Z_c'	858	'L'
'tag1601819Z1835Z9Z12Z_c'	854	'L'
'tag1601819Z1845Z9Z12Z15Z_c'	852	'L'
'tag1801819Z160_c'	860	'L'
'tag1801819Z180_c'	888	'L'
'tag1801819Z18111Z_c'	886	'L'
'tag1801819Z1819Z_c'	886	'L'
'tag1801819Z1835Z9Z12Z_c'	882	'L'
'tag1801819Z1845Z9Z12Z15Z_c'	880	'L'
'tag18111Z18111Z160_c'	858	'L'
'tag18111Z18111Z180_c'	886	'L'
'tag18111Z18111Z18111Z_c'	884	'L'
'tag18111Z18111Z1819Z_c'	884	'L'
'tag18111Z18111Z1835Z9Z12Z_c'	880	'L'
'tag18111Z18111Z1845Z9Z12Z15Z_c'	878	'L'
'tag18111Z1819Z160_c'	858	'L'
'tag18111Z1819Z180_c'	886	'L'
'tag18111Z1819Z18111Z_c'	884	'L'
'tag18111Z1819Z1819Z_c'	884	'L'
'tag18111Z1819Z1835Z9Z12Z_c'	880	'L'
'tag18111Z1819Z1845Z9Z12Z15Z_c'	878	'L'
'tag1819Z18111Z160_c'	858	'L'
'tag1819Z18111Z180_c'	886	'L'
'tag1819Z18111Z18111Z_c'	884	'L'
'tag1819Z18111Z1819Z_c'	884	'L'
'tag1819Z18111Z1835Z9Z12Z_c'	880	'L'
'tag1819Z18111Z1845Z9Z12Z15Z_c'	878	'L'
'tag1819Z1819Z160_c'	858	'L'
'tag1819Z1819Z180_c'	886	'L'
'tag1819Z1819Z18111Z_c'	884	'L'
'tag1819Z1819Z1819Z_c'	884	'L'
'tag1819Z1819Z1835Z9Z12Z_c'	880	'L'
'tag1819Z1819Z1845Z9Z12Z15Z_c'	878	'L'
'thmmp_c'	343	'N'
'thrtrna_c'	429	'P'
'trnaala_c'	328	'P'
'trnaarg_c'	462	'P'
'trnaasn_c'	328	'P'
'trnaasp_c'	328	'P'
'trnacys_c'	462	'P'
'trnagln_c'	462	'P'
'trnaglu_c'	462	'P'
'trnagly_c'	328	'P'
'trnahis_c'	328	'P'
'trnaile_c'	462	'P'
'trnaleu_c'	462	'P'
'trnalys_c'	328	'P'
'trnamet_c'	462	'P'
'trnaphe_c'	328	'P'
'trnapro_c'	328	'P'
'trnaser_c'	328	'P'
'trnathr_c'	328	'P'
'trnatrp_c'	462	'P'
'trnatyr_c'	462	'P'
'trnaval_c'	462	'P'
'trptrna_c'	648	'P'
'tyrtrna_c'	625	'P'
'utp_c'	480	'N'
'valtrna_c'	561	'P'
'vioxan_u'	600	'N'
'zaxan_u'	568	'N'};




%Get main fractions:
[P,X] = getFraction(model,comps,'P',0);
[C,X] = getFraction(model,comps,'C',X);
[R,X] = getFraction(model,comps,'R',X);
[D,X] = getFraction(model,comps,'D',X);
[L,X] = getFraction(model,comps,'L',X);

%Add up any remaining components:
bioPos = strcmp(model.rxns,'r_4041');
for i = 1:length(model.mets)
    pos = strcmp(comps(:,1),model.mets{i});
    if sum(pos) == 1
        abundance = -model.S(i,bioPos)*comps{pos,2}/1000;
        X         = X + abundance;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,X] = getFraction(model,comps,compType,X)

%Define pseudoreaction name:
rxnName = [compType ' pseudoreaction'];
rxnName = strrep(rxnName,'P','protein');
rxnName = strrep(rxnName,'C','carbohydrate');
rxnName = strrep(rxnName,'R','RNA');
rxnName = strrep(rxnName,'D','DNA');
rxnName = strrep(rxnName,'L','lipid backbone');

%Add up fraction:
fractionPos = strcmp(model.rxnNames,rxnName);
if contains(rxnName,'lipid')
    subs = model.S(:,fractionPos) < 0;        %substrates in pseudo-rxn
    F    = -sum(model.S(subs,fractionPos));   %g/gDW
else
    comps = comps(strcmp(comps(:,3),compType),:);
    F = 0;
    %Add up all components:
    for i = 1:length(model.mets)
        pos = strcmp(comps(:,1),model.mets{i});
        if sum(pos) == 1
            abundance = -model.S(i,fractionPos)*(comps{pos,2}-18)/1000;
            F         = F + abundance;
        end
    end
end
X = X + F;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
