const int nBTagAnalyzerTaggers = 37;
TString BTagAnalyzerTaggerName[nBTagAnalyzerTaggers] = {"TCHPL", "TCHET", "TCHPM", "TCHPT", "TCHEL", "SSVHEL", 
							"CSVL",  "CSVM",  "CSVT", 
							"JBPT",  "JBPT1", "JBPM",  "JBPM1", "JBPM2", "JBPM3", "JBPM4", 
							"JBPM5", "JBPM6", "JBPM7", "JBPM8", "JBPM9", "JBPL", 
							"JPL",   "JPM",   "JPT", 
							"CSVv2L","CSVv2M","CSVv2T", 
							"CMVAv2L", "CMVAv2M", "CMVAv2T",
							"DeepCSVL", "DeepCSVM", "DeepCSVT",
							"DeepFlavourL", "DeepFlavourM", "DeepFlavourT"};
float   BTagAnalyzerTaggerCut[nBTagAnalyzerTaggers]  = {1.19,    10.2,    1.93,    3.41,    1.7,     3.05,
							0.244,   0.679,   0.898,
							3.74,    3.14,    2.55,    2.43,    2.31,    2.18,    2.06,
							1.94,    1.82,    1.70,    1.57,    1.45,    1.33,
							0.2620,  0.5344,  0.7914,
							0.5803,  0.8838,  0.9693,
						       -0.5884,  0.4432,  0.9432,
						        0.2217,  0.6321,  0.8953,
						        0.0614,  0.3093,  0.7221};
