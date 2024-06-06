############Figure 4A############

library(multtest)
library(Seurat)
library(dplyr)
library(mindr)
library(tidyverse)
DimPlot(Tcell, reduction = "umap", label = TRUE,cols = mycolor, raster=FALSE) + NoLegend()

############Figure 4B############

DimPlot(Tcell, cells.highlight = Cells(subset(Tcell, 
                                              subset = cdr3 == "CASSDSGDTNTGELFF"|cdr3 == "CASSELAGGLETQYF"|
                                                cdr3 == "CAKTLAKNIQYF"|cdr3 == "CASSKGETDIQYF"|cdr3 == "CASSVGDTDTQYF"|
                                                cdr3 == "CASSGEGGTSNEQFF"|cdr3 == "CATSDSTGVETQYF"|
                                                cdr3 == "CASSLQGENGELFF"|cdr3 == "CSVVGREGQPQHF"|cdr3 == "CASSIGRGDTQYF"|
                                                cdr3 == "CASSPGLQETQYF"|cdr3 == "CASSRTDGGRGYTF"|
                                                cdr3 == "CASSPGQSSGANVLTF"|cdr3 == "CASSIDRGSDTQYF"|
                                                cdr3 == "CASKAGTGDDTGELFF"|cdr3 == "CASSSGSGGTDTQYF"|cdr3 == "CASSPGTSSYNEQFF"|
                                                cdr3 == "CASSIPGADYNEQFF"|cdr3 == "CASSQRPPGEDSYEQYF"|cdr3 == "CASSVQGASPQHF"|
                                                cdr3 == "CASKTGGRTGELFF"|cdr3 == "CASSERAGEADTQYF"|cdr3 == "CASSEQRATDTQYF"|
                                                cdr3 == "CASSDAAGGGTDTQYF"|cdr3 == "CASSIRVSSGNTIYF"|cdr3 == "CAISGGEEGHEQYF"|
                                                cdr3 == "CSAREPADYNSPLHF"|cdr3 == "CATTFREGQDYNEQFF"|cdr3 == "CASSIPGTDTQYF"|
                                                cdr3 == "CASSPSRGGDEKLFF"|cdr3 == "CASSLGIRSDQPQHF"|
                                                cdr3 == "CASRVSNQPQHF"|cdr3 == "CASRRDSYEQFF"|cdr3 == "CASSLAEPHEQYV"|
                                                cdr3 == "CSARDPGTGRNTEAFF"|cdr3 == "CSVPGQGDTGELFF"|cdr3 == "CASRVNNVLTF"|
                                                cdr3 == "CASSQGDLAGTQYF"|cdr3 == "CASSPPGQGASEQFF"|cdr3 == "CASSRQGQGTEAFF"|
                                                cdr3 == "CASSSLPGGGNSNQPQHF"|cdr3 == "CASSLIPGGGDSPLHF"|cdr3 == "CASSRYPSGPAETQYF"|
                                                cdr3 == "CASSESSGDYNEQFF"|cdr3 == "CASSFRQTLNTEAFF"|
                                                cdr3 == "CAISETEETQYF"|cdr3 == "CSARDWQLGNRGYGYTF"|cdr3 == "CSARSPDRGTEAFF"|
                                                cdr3 == "CSALQGDLSYNEQFF"|cdr3 == "CSAVTSGGTDTQYF"|cdr3 == "CASSQEGTGSGSISYNEQFF"|
                                                cdr3 == "CASSQALRDYNTGELFF"|cdr3 == "CASSPTDSPVETQYF"|cdr3 == "CASSLYTGDQPQHF"|
                                                cdr3 == "CASSHGPTGYTF"|cdr3 == "CASIEGKTEAFF"|cdr3 == "CATHEASVGNEQFF"|
                                                cdr3 == "CASSQDGDRAGQPQHF"|cdr3 == "CASSRQGMGGNQPQHF"|cdr3 == "CASSLTTGTSGTHTDTQYF"|
                                                cdr3 == "CASSLGAGSAYEQYF"|cdr3 == "CASSGGQYTEAFF"|cdr3 == "CASSSNLRDTGKNTEAFF"|
                                                cdr3 == "CASSLAGGNQPQHF"|cdr3 == "CASSFRDRVTNEQFF"|cdr3 == "CASRGTDSGNTIYF"|
                                                cdr3 == "CASRGDEDEQYV"|cdr3 == "CASSDSPSGGGTDTQYF"|cdr3 == "CASSGTETETQYF"|
                                                cdr3 == "CASSYPPAGGADEQFF"|cdr3 == "CASSYGQENNQPQHF"|cdr3 == "CASSYGGAGANVLTF"|
                                                cdr3 == "CASSSGQTGWIDTDTQYF"|cdr3 == "CASSWGRDQETQYF"|cdr3 == "CASGDRGRVLRDSPLHF"|
                                                cdr3 == "CASGGTGTQPQHF"|cdr3 == "CASCTGTPTGELFF"|cdr3 == "CASQSGAYNEQFF"|
                                                
                                                cdr3 =="CASSLQTGTIGGYTF"|cdr3 =="CASSQDRGGGSTDTQYF"|cdr3 =="CASSEGGGSYNEQFF"|
                                                cdr3 =="CAISEGPFSRGHRNEQFF"|cdr3 =="CASSQGQRYEQYF"|
                                                cdr3 =="CSVENPLDRGDYGYTF"|cdr3 =="CASSAPGTGVEQYF"|
                                                cdr3 =="CASDRGGYEQYF"|cdr3 =="CASSPYTATNEKLFF"|cdr3 =="CASTQNTGELFF"|
                                                cdr3 =="CASSDSGGGAGDTQYF"|cdr3 =="CSVEDPDQVGGELFF"|cdr3 =="CASSSGPLGGYTF"|
                                                cdr3 =="CASMDRVGGTDTQYF"|cdr3 =="CASSQDYPRGEQYF"|cdr3 =="CASSEASGGAETQYF"|
                                                cdr3 =="CASSSPDRESGELFF"|cdr3 =="CASSLGTSGGSETQYF"|cdr3 =="CASSLGGPEQYV"|
                                                cdr3 =="CASSPGGGQPQHF"|cdr3 =="CSVEEGRAVGEAFF"|
                                                cdr3 =="CATSGTSGSYNEQFF"|cdr3 =="CASSRQGHYEQYF"|cdr3 =="CASSSRTGEETQYF"|
                                                cdr3 =="CSARDRNRGPNYGYTF"|cdr3 =="CASSQDGYRGNYGYTF"|cdr3 =="CATSSDDGGLEEKLFF"|
                                                cdr3 =="CASSQVQLRDRANYGYTF"|cdr3 =="CASSSFSTQQYF"|cdr3 =="CASSLTQGQETQYF"|
                                                cdr3 =="CASSIQGRRETQYF"|cdr3 =="CASSEAVNLGEQYF"|cdr3 =="CASSYKQEQPQHF"|
                                                cdr3 =="CASSLTGIYEQYF"|cdr3 =="CASSFGEKGPYEQYF"|cdr3 =="CASSFKGDRGGVAFF"|
                                                cdr3 =="CASAGFEQFF"|cdr3 =="CSARVGHLDTGELFF"|cdr3 =="CASKTGGGDQPQHF"|
                                                cdr3 =="CASTGGGSYEQYV"|cdr3 =="CASSQVGGSTDTQYF"|cdr3 =="CASSPGNTEAFF"|
                                                cdr3 =="CASSDRANYGYTF"|cdr3 =="CASSLRGAGANVLTF"|cdr3 =="CASSLRGGSSYEQYF"|
                                                cdr3 =="CASSIRTLEAQHF"|cdr3 =="CASSYSHRRYSGANVLTF"|cdr3 =="CASSSTSGGATDTQYF"|
                                                cdr3 =="CASSLGGSLNTQYF"|
                                                
                                                cdr3 == "CASGTTDTQYF"|cdr3 == "CASSLESGGEQFF"|
                                                cdr3 == "CASSLEGRLGDTQYF"|cdr3 == "CASSTDRNSPLHF"|cdr3 == "CASTNRGNQPQHF"|
                                                cdr3 == "CASSRKSSYNEQFF"|cdr3 == "CASSLSPGSNSPLHF"|cdr3 == "CSVVNPLNSPYEQYF"|
                                                cdr3 == "CASRRLAGQQETQYF"|cdr3 == "CASSVDGGTDTQYF"|cdr3 == "CAARTTGPGNTIYF"|
                                                cdr3 == "CSARAFGRFSTDTQYF"|cdr3 == "CASSLDGGGTEAFF"|
                                                cdr3 == "CASSLEPSGSGEQYF"|cdr3 == "CASSQADRGDEQYF"|cdr3 == "CASSHGARYNEQFF"|
                                                cdr3 == "CASSPADQETQYF"|cdr3 == "CASSLAPLGADYNEQFF"|cdr3 == "CASSPLAGGNTQYF"|
                                                cdr3 == "CASSLSGGPGETQYF"|cdr3 == "CASSLDGMNTEAFF"|cdr3 == "CASSEMRNNEQFF"|
                                                cdr3 == "CASSFSGRAGSEQFF"|cdr3 == "CSVGGNSVETQYF"|
                                                cdr3 == "CSAPWTGGRYEQYF"|cdr3 == "CASSLEGTSNNEQFF"|cdr3 == "CATSRDRSSGSEQYF"|
                                                cdr3 == "CASSQEAALRTAQYF"|cdr3 == "CASSLGRQSGDRAFF"|cdr3 == "CASSNQGRNQPQHF"|
                                                cdr3 == "CASSPGVAGGTGELFF"|cdr3 == "CASSGSGNQPQHF"|cdr3 == "CASRGRTEGSYNEQFF"|
                                                cdr3 == "CASSPGFSKETQYF"|cdr3 == "CASSLEPAGVEETQYF"|cdr3 == "CASSRAGGSSYEQYF"|
                                                cdr3 == "CATRASGKMDEQFF"|cdr3 == "CATSRAWSSYNEQFF"|cdr3 == "CASTPRRLAKNIQYF"|
                                                cdr3 == "CASSPPGHTEAFF"|cdr3 == "CASSETGALGAGASHEQYF"|cdr3 == "CASSLEGGGASGYTF"|
                                                cdr3 == "CASSEYGGAGQETQYF"|cdr3 == "CASSSSPNQPQHF"|cdr3 == "CASSWDRGNEQFF"|
                                                cdr3 == "CASSPQGGTEAFF"|
                                                
                                                cdr3 =="CASSQSGTSGDNEQFF"|cdr3 =="CASTLDYNEQFF"|cdr3 =="CASSEGQGYNQPQHF"|
                                                cdr3 =="CASSFWARGQGFASTDTQYF"|cdr3 =="CASSRDPRTSGPLPGDTQYF"|
                                                cdr3 =="CSARGGRGDTEAFF"|cdr3 =="CASSTGTSGYNEQYF"|cdr3 =="CASSLRDSSNQPQHF"|
                                                cdr3 =="CASSPNPRTSGKDEQFF"|cdr3 =="CASSQDRAADGYTF"|cdr3 =="CASSLWGAETQYF"|
                                                cdr3 =="CASSLGTSNEQFF"|cdr3 =="CATSRFAGEEQYF"|cdr3 =="CASSHASGNFYEQYF"|
                                                cdr3 =="CASSPPQGSNTEAFF"|cdr3 =="CASSLEGGLSSEAFF"|cdr3 =="CASSFVAGGRNEQFF"|
                                                cdr3 =="CASQGNTDTQYF"|cdr3 =="CASTSLTSEAFF"|cdr3 =="CASTLAGGSTYEQYF"|
                                                cdr3 =="CASSLARDNEQFF"|cdr3 =="CASSLGPWNEQFF"|cdr3 =="CSADSGAGTGELFF"|
                                                cdr3 =="CSARDVAGTSYNEQFF"|cdr3 =="CSAAREGNTEAFF"|cdr3 =="CASSLAGNNYGYTF"|
                                                cdr3 =="CASSDSTSGSTDTQYF"|cdr3 =="CASSLTSDRGLGNTIYF"|cdr3 =="CASSLASWAAEQYF"|
                                                cdr3 =="CASRGVRAYGYTF"|cdr3 =="CASSEPRTNYNEQFF"|cdr3 =="CASSEDRVGEQYF"|
                                                cdr3 =="CAISIDFHEQYF"|cdr3 =="CAISGLTNSYNEQFF"|cdr3 =="CASSGATRVSYEQYF"|
                                                cdr3 =="CASSKSGGNNEQFF"|cdr3 =="CSVRRGGNEQFF"|cdr3 =="CSVEVGTSGSTEQFF"|
                                                cdr3 =="CSAGRASSYEQYF"|cdr3 =="CSARESGRRGRPYEQYF"|cdr3 =="CSARDGRLFYEQYF"|
                                                cdr3 =="CASSLEQIVRTDTQYF"|cdr3 =="CASSLEQGVRSEQFF"|cdr3 =="CASSPWGGQGTDTQYF"|
                                                cdr3 =="CASSLKGPAYNEQFF"|cdr3 =="CASSLDLGRSYEQYF"|cdr3 =="CASSRRENSPLHF"|
                                                cdr3 =="CASRESADHNEQFF"|cdr3 =="CSVRQGGNEQFF"|cdr3 =="CSAFPSRGAGQPQHF"|
                                                cdr3 =="CSARDLASGYNEQFF"|cdr3 =="CATSRAAGGPPDTQYF"|cdr3 =="CASSQSGVGRGYEQYF"|
                                                cdr3 =="CASSPAWIAGEDTQYF"|cdr3 =="CASRPQGYSNQPQHF"|
                                                cdr3 =="CASSKTGPTKTYEQYF"|cdr3 =="CASSYVDATEETQYF"|cdr3 =="CASRRWTGDDEQYF"|
                                                
                                                cdr3 =="CASSRHSGLFLYEQYF"|cdr3 =="CASSLGPEASNEQFF"|cdr3 =="CASSLSGGEGPNEQFF"|
                                                cdr3 =="CASQPRADNSPLHF"|cdr3 =="CASSSLVQETQYF"|cdr3 =="CASSPLSGEQYF"|
                                                cdr3 =="CASSVEGGTDTGELFF"|cdr3 =="CASSSRGDNSPLHF"|cdr3 =="CASSPRDRNIQYF"|
                                                cdr3 =="CASSPVTGTGNYGYTF"|cdr3 =="CASSVALNLAGVYNEQFF"|cdr3 =="CASSPGQGRNEKLFF"|
                                                cdr3 =="CSASLGPGELFF"|cdr3 =="CATGGGYTF"|cdr3 =="CASSLGLAGVNEQYF"|
                                                cdr3 =="CSARTEKALSIQYF"|cdr3 =="CASSVSVGEQYF"|cdr3 =="CSASTTSGTPSTQYF"|
                                                cdr3 =="CAGGGGDGNIQYF"|cdr3 =="CASSVFGGTGGEETQYF"|cdr3 =="CATEPGIGELFF"|
                                                cdr3 =="CASSLAGTDTQYF"|cdr3 =="CASGLSTLAETQYF"|cdr3 =="CSAPTLLGEQFF"|
                                                cdr3 =="CASSLRVPSSTDTQYF"|cdr3 =="CASKLTGTSGNEQFF"|cdr3 =="CASTLEGSLNTEAFF"|
                                                cdr3 =="CASSQQGVMSNQPQHF"|cdr3 =="CASSVPENLNNEQFF"|cdr3 =="CSVVGPWGEGSIYEQYF"|
                                                cdr3 =="CSVDSPSRSTDTQYF"|cdr3 =="CASRDNRGSGNTIYF")),
        split.by = "sample",raster=FALSE,pt.size = 1.5)


############Figure 4C-D############
TCR_react <- subset(Tcell, 
                    subset = cdr3 == "CASSDSGDTNTGELFF"|cdr3 == "CASSELAGGLETQYF"|
                      cdr3 == "CAKTLAKNIQYF"|cdr3 == "CASSKGETDIQYF"|cdr3 == "CASSVGDTDTQYF"|
                      cdr3 == "CASSGEGGTSNEQFF"|cdr3 == "CATSDSTGVETQYF"|
                      cdr3 == "CASSLQGENGELFF"|cdr3 == "CSVVGREGQPQHF"|cdr3 == "CASSIGRGDTQYF"|
                      cdr3 == "CASSPGLQETQYF"|cdr3 == "CASSRTDGGRGYTF"|
                      cdr3 == "CASSPGQSSGANVLTF"|cdr3 == "CASSIDRGSDTQYF"|
                      cdr3 == "CASKAGTGDDTGELFF"|cdr3 == "CASSSGSGGTDTQYF"|cdr3 == "CASSPGTSSYNEQFF"|
                      cdr3 == "CASSIPGADYNEQFF"|cdr3 == "CASSQRPPGEDSYEQYF"|cdr3 == "CASSVQGASPQHF"|
                      cdr3 == "CASKTGGRTGELFF"|cdr3 == "CASSERAGEADTQYF"|cdr3 == "CASSEQRATDTQYF"|
                      cdr3 == "CASSDAAGGGTDTQYF"|cdr3 == "CASSIRVSSGNTIYF"|cdr3 == "CAISGGEEGHEQYF"|
                      cdr3 == "CSAREPADYNSPLHF"|cdr3 == "CATTFREGQDYNEQFF"|cdr3 == "CASSIPGTDTQYF"|
                      cdr3 == "CASSPSRGGDEKLFF"|cdr3 == "CASSLGIRSDQPQHF"|
                      cdr3 == "CASRVSNQPQHF"|cdr3 == "CASRRDSYEQFF"|cdr3 == "CASSLAEPHEQYV"|
                      cdr3 == "CSARDPGTGRNTEAFF"|cdr3 == "CSVPGQGDTGELFF"|cdr3 == "CASRVNNVLTF"|
                      cdr3 == "CASSQGDLAGTQYF"|cdr3 == "CASSPPGQGASEQFF"|cdr3 == "CASSRQGQGTEAFF"|
                      cdr3 == "CASSSLPGGGNSNQPQHF"|cdr3 == "CASSLIPGGGDSPLHF"|cdr3 == "CASSRYPSGPAETQYF"|
                      cdr3 == "CASSESSGDYNEQFF"|cdr3 == "CASSFRQTLNTEAFF"|
                      cdr3 == "CAISETEETQYF"|cdr3 == "CSARDWQLGNRGYGYTF"|cdr3 == "CSARSPDRGTEAFF"|
                      cdr3 == "CSALQGDLSYNEQFF"|cdr3 == "CSAVTSGGTDTQYF"|cdr3 == "CASSQEGTGSGSISYNEQFF"|
                      cdr3 == "CASSQALRDYNTGELFF"|cdr3 == "CASSPTDSPVETQYF"|cdr3 == "CASSLYTGDQPQHF"|
                      cdr3 == "CASSHGPTGYTF"|cdr3 == "CASIEGKTEAFF"|cdr3 == "CATHEASVGNEQFF"|
                      cdr3 == "CASSQDGDRAGQPQHF"|cdr3 == "CASSRQGMGGNQPQHF"|cdr3 == "CASSLTTGTSGTHTDTQYF"|
                      cdr3 == "CASSLGAGSAYEQYF"|cdr3 == "CASSGGQYTEAFF"|cdr3 == "CASSSNLRDTGKNTEAFF"|
                      cdr3 == "CASSLAGGNQPQHF"|cdr3 == "CASSFRDRVTNEQFF"|cdr3 == "CASRGTDSGNTIYF"|
                      cdr3 == "CASRGDEDEQYV"|cdr3 == "CASSDSPSGGGTDTQYF"|cdr3 == "CASSGTETETQYF"|
                      cdr3 == "CASSYPPAGGADEQFF"|cdr3 == "CASSYGQENNQPQHF"|cdr3 == "CASSYGGAGANVLTF"|
                      cdr3 == "CASSSGQTGWIDTDTQYF"|cdr3 == "CASSWGRDQETQYF"|cdr3 == "CASGDRGRVLRDSPLHF"|
                      cdr3 == "CASGGTGTQPQHF"|cdr3 == "CASCTGTPTGELFF"|cdr3 == "CASQSGAYNEQFF"|
                      
                      cdr3 =="CASSLQTGTIGGYTF"|cdr3 =="CASSQDRGGGSTDTQYF"|cdr3 =="CASSEGGGSYNEQFF"|
                      cdr3 =="CAISEGPFSRGHRNEQFF"|cdr3 =="CASSQGQRYEQYF"|
                      cdr3 =="CSVENPLDRGDYGYTF"|cdr3 =="CASSAPGTGVEQYF"|
                      cdr3 =="CASDRGGYEQYF"|cdr3 =="CASSPYTATNEKLFF"|cdr3 =="CASTQNTGELFF"|
                      cdr3 =="CASSDSGGGAGDTQYF"|cdr3 =="CSVEDPDQVGGELFF"|cdr3 =="CASSSGPLGGYTF"|
                      cdr3 =="CASMDRVGGTDTQYF"|cdr3 =="CASSQDYPRGEQYF"|cdr3 =="CASSEASGGAETQYF"|
                      cdr3 =="CASSSPDRESGELFF"|cdr3 =="CASSLGTSGGSETQYF"|cdr3 =="CASSLGGPEQYV"|
                      cdr3 =="CASSPGGGQPQHF"|cdr3 =="CSVEEGRAVGEAFF"|
                      cdr3 =="CATSGTSGSYNEQFF"|cdr3 =="CASSRQGHYEQYF"|cdr3 =="CASSSRTGEETQYF"|
                      cdr3 =="CSARDRNRGPNYGYTF"|cdr3 =="CASSQDGYRGNYGYTF"|cdr3 =="CATSSDDGGLEEKLFF"|
                      cdr3 =="CASSQVQLRDRANYGYTF"|cdr3 =="CASSSFSTQQYF"|cdr3 =="CASSLTQGQETQYF"|
                      cdr3 =="CASSIQGRRETQYF"|cdr3 =="CASSEAVNLGEQYF"|cdr3 =="CASSYKQEQPQHF"|
                      cdr3 =="CASSLTGIYEQYF"|cdr3 =="CASSFGEKGPYEQYF"|cdr3 =="CASSFKGDRGGVAFF"|
                      cdr3 =="CASAGFEQFF"|cdr3 =="CSARVGHLDTGELFF"|cdr3 =="CASKTGGGDQPQHF"|
                      cdr3 =="CASTGGGSYEQYV"|cdr3 =="CASSQVGGSTDTQYF"|cdr3 =="CASSPGNTEAFF"|
                      cdr3 =="CASSDRANYGYTF"|cdr3 =="CASSLRGAGANVLTF"|cdr3 =="CASSLRGGSSYEQYF"|
                      cdr3 =="CASSIRTLEAQHF"|cdr3 =="CASSYSHRRYSGANVLTF"|cdr3 =="CASSSTSGGATDTQYF"|
                      cdr3 =="CASSLGGSLNTQYF"|
                      
                      cdr3 == "CASGTTDTQYF"|cdr3 == "CASSLESGGEQFF"|
                      cdr3 == "CASSLEGRLGDTQYF"|cdr3 == "CASSTDRNSPLHF"|cdr3 == "CASTNRGNQPQHF"|
                      cdr3 == "CASSRKSSYNEQFF"|cdr3 == "CASSLSPGSNSPLHF"|cdr3 == "CSVVNPLNSPYEQYF"|
                      cdr3 == "CASRRLAGQQETQYF"|cdr3 == "CASSVDGGTDTQYF"|cdr3 == "CAARTTGPGNTIYF"|
                      cdr3 == "CSARAFGRFSTDTQYF"|cdr3 == "CASSLDGGGTEAFF"|
                      cdr3 == "CASSLEPSGSGEQYF"|cdr3 == "CASSQADRGDEQYF"|cdr3 == "CASSHGARYNEQFF"|
                      cdr3 == "CASSPADQETQYF"|cdr3 == "CASSLAPLGADYNEQFF"|cdr3 == "CASSPLAGGNTQYF"|
                      cdr3 == "CASSLSGGPGETQYF"|cdr3 == "CASSLDGMNTEAFF"|cdr3 == "CASSEMRNNEQFF"|
                      cdr3 == "CASSFSGRAGSEQFF"|cdr3 == "CSVGGNSVETQYF"|
                      cdr3 == "CSAPWTGGRYEQYF"|cdr3 == "CASSLEGTSNNEQFF"|cdr3 == "CATSRDRSSGSEQYF"|
                      cdr3 == "CASSQEAALRTAQYF"|cdr3 == "CASSLGRQSGDRAFF"|cdr3 == "CASSNQGRNQPQHF"|
                      cdr3 == "CASSPGVAGGTGELFF"|cdr3 == "CASSGSGNQPQHF"|cdr3 == "CASRGRTEGSYNEQFF"|
                      cdr3 == "CASSPGFSKETQYF"|cdr3 == "CASSLEPAGVEETQYF"|cdr3 == "CASSRAGGSSYEQYF"|
                      cdr3 == "CATRASGKMDEQFF"|cdr3 == "CATSRAWSSYNEQFF"|cdr3 == "CASTPRRLAKNIQYF"|
                      cdr3 == "CASSPPGHTEAFF"|cdr3 == "CASSETGALGAGASHEQYF"|cdr3 == "CASSLEGGGASGYTF"|
                      cdr3 == "CASSEYGGAGQETQYF"|cdr3 == "CASSSSPNQPQHF"|cdr3 == "CASSWDRGNEQFF"|
                      cdr3 == "CASSPQGGTEAFF"|
                      
                      cdr3 =="CASSQSGTSGDNEQFF"|cdr3 =="CASTLDYNEQFF"|cdr3 =="CASSEGQGYNQPQHF"|
                      cdr3 =="CASSFWARGQGFASTDTQYF"|cdr3 =="CASSRDPRTSGPLPGDTQYF"|
                      cdr3 =="CSARGGRGDTEAFF"|cdr3 =="CASSTGTSGYNEQYF"|cdr3 =="CASSLRDSSNQPQHF"|
                      cdr3 =="CASSPNPRTSGKDEQFF"|cdr3 =="CASSQDRAADGYTF"|cdr3 =="CASSLWGAETQYF"|
                      cdr3 =="CASSLGTSNEQFF"|cdr3 =="CATSRFAGEEQYF"|cdr3 =="CASSHASGNFYEQYF"|
                      cdr3 =="CASSPPQGSNTEAFF"|cdr3 =="CASSLEGGLSSEAFF"|cdr3 =="CASSFVAGGRNEQFF"|
                      cdr3 =="CASQGNTDTQYF"|cdr3 =="CASTSLTSEAFF"|cdr3 =="CASTLAGGSTYEQYF"|
                      cdr3 =="CASSLARDNEQFF"|cdr3 =="CASSLGPWNEQFF"|cdr3 =="CSADSGAGTGELFF"|
                      cdr3 =="CSARDVAGTSYNEQFF"|cdr3 =="CSAAREGNTEAFF"|cdr3 =="CASSLAGNNYGYTF"|
                      cdr3 =="CASSDSTSGSTDTQYF"|cdr3 =="CASSLTSDRGLGNTIYF"|cdr3 =="CASSLASWAAEQYF"|
                      cdr3 =="CASRGVRAYGYTF"|cdr3 =="CASSEPRTNYNEQFF"|cdr3 =="CASSEDRVGEQYF"|
                      cdr3 =="CAISIDFHEQYF"|cdr3 =="CAISGLTNSYNEQFF"|cdr3 =="CASSGATRVSYEQYF"|
                      cdr3 =="CASSKSGGNNEQFF"|cdr3 =="CSVRRGGNEQFF"|cdr3 =="CSVEVGTSGSTEQFF"|
                      cdr3 =="CSAGRASSYEQYF"|cdr3 =="CSARESGRRGRPYEQYF"|cdr3 =="CSARDGRLFYEQYF"|
                      cdr3 =="CASSLEQIVRTDTQYF"|cdr3 =="CASSLEQGVRSEQFF"|cdr3 =="CASSPWGGQGTDTQYF"|
                      cdr3 =="CASSLKGPAYNEQFF"|cdr3 =="CASSLDLGRSYEQYF"|cdr3 =="CASSRRENSPLHF"|
                      cdr3 =="CASRESADHNEQFF"|cdr3 =="CSVRQGGNEQFF"|cdr3 =="CSAFPSRGAGQPQHF"|
                      cdr3 =="CSARDLASGYNEQFF"|cdr3 =="CATSRAAGGPPDTQYF"|cdr3 =="CASSQSGVGRGYEQYF"|
                      cdr3 =="CASSPAWIAGEDTQYF"|cdr3 =="CASRPQGYSNQPQHF"|
                      cdr3 =="CASSKTGPTKTYEQYF"|cdr3 =="CASSYVDATEETQYF"|cdr3 =="CASRRWTGDDEQYF"|
                      
                      cdr3 =="CASSRHSGLFLYEQYF"|cdr3 =="CASSLGPEASNEQFF"|cdr3 =="CASSLSGGEGPNEQFF"|
                      cdr3 =="CASQPRADNSPLHF"|cdr3 =="CASSSLVQETQYF"|cdr3 =="CASSPLSGEQYF"|
                      cdr3 =="CASSVEGGTDTGELFF"|cdr3 =="CASSSRGDNSPLHF"|cdr3 =="CASSPRDRNIQYF"|
                      cdr3 =="CASSPVTGTGNYGYTF"|cdr3 =="CASSVALNLAGVYNEQFF"|cdr3 =="CASSPGQGRNEKLFF"|
                      cdr3 =="CSASLGPGELFF"|cdr3 =="CATGGGYTF"|cdr3 =="CASSLGLAGVNEQYF"|
                      cdr3 =="CSARTEKALSIQYF"|cdr3 =="CASSVSVGEQYF"|cdr3 =="CSASTTSGTPSTQYF"|
                      cdr3 =="CAGGGGDGNIQYF"|cdr3 =="CASSVFGGTGGEETQYF"|cdr3 =="CATEPGIGELFF"|
                      cdr3 =="CASSLAGTDTQYF"|cdr3 =="CASGLSTLAETQYF"|cdr3 =="CSAPTLLGEQFF"|
                      cdr3 =="CASSLRVPSSTDTQYF"|cdr3 =="CASKLTGTSGNEQFF"|cdr3 =="CASTLEGSLNTEAFF"|
                      cdr3 =="CASSQQGVMSNQPQHF"|cdr3 =="CASSVPENLNNEQFF"|cdr3 =="CSVVGPWGEGSIYEQYF"|
                      cdr3 =="CSVDSPSRSTDTQYF"|cdr3 =="CASRDNRGSGNTIYF" )
table(Idents(TCR_react), TCR_react$sample)#
TCR_react_percent = prop.table(table(TCR_react$sample,Idents(TCR_react)), margin = 1)
TCR_react_count = table(TCR_react$sample,Idents(TCR_react))
TCR_react_percent
write.csv(TCR_react_percent,"TCR_react_percent_sample.csv")
write.csv(TCR_react_count, "TCR_react_count_sample.csv")

############Figure 4E-F############
library(ggsci)
mycolor <- c(pal_npg()(10),pal_jama()(7), pal_uchicago()(7), pal_jama()(7))

library(RColorBrewer)
library(slingshot)
library(SingleCellExperiment)

pbmc.sce <- as.SingleCellExperiment(CD8)
sim <- slingshot(pbmc.sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP',   start.clus = 1 )
plot(reducedDims(sim)$UMAP, col = mycolor[sim$seurat_clusters], pch=16, asp = 2)
lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black')

pbmc.sce <- as.SingleCellExperiment(CD4)
sim <- slingshot(pbmc.sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP',   start.clus = 1 )
plot(reducedDims(sim)$UMAP, col = mycolor[sim$seurat_clusters], pch=16, asp = 2)
lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black')

############Figure 4G############
library(msigdbr) 
library(fgsea)
library(Seurat)

DefaultAssay(Tcell) <- "RNA"
markers <- FindMarkers(Tcell, ident.1 = "CD8_Cytotoxic", min.pct = 0.25, logfc.threshold = 0)  
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes, avg_log2FC)
ranks<- deframe(cluster.genes)
####GSEA####
mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
fgsea_sets = mdb_c2 [grep("KEGG",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
write.csv(fgseaRes[,1:7], "CD8_GZMB_KEGG.csv")

markers <- FindMarkers(Tcell, ident.1 = "CD4_Cytotoxic", min.pct = 0.25, logfc.threshold = 0)  
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes, avg_log2FC)
ranks<- deframe(cluster.genes)
####GSEA####
mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
fgsea_sets = mdb_c2 [grep("KEGG",mdb_c2 $gs_name),] %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
write.csv(fgseaRes[,1:7], "CD4_GZMB_KEGG.csv")


############Figure 4H############
markers = c("CX3CR1","CXCR1","CXCR2","CXCR3","CXCR4","CXCR5","CXCR6","CCR1","CCR2"
            ,"CCR3","CCR4", "CCR5","CCR6",'CCR7', 'CCR8',"CCR9", "CCR10","XCR1")
DotPlot(Tcell, features = markers,
        assay='RNA' ) + 
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 


############Figure 4I############
VlnPlot(pbmc, features = "CX3CL1")
