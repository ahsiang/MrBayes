#ifndef __PROPOSAL_H__
#define __PROPOSAL_H__

int     Move_Aamodel (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_AddBranch (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Adgamma (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_AddDeleteCPPEvent (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Allocation (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Alphadir_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Beta (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Beta_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_BrLen (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ClockRate_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_CPPEventPosition (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_CPPRate (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_CPPRateMultiplier_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_CPPRateMultiplierRnd (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_DelBranch (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Extinction (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ExtFossilSPRClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ExtSPR (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ExtSPR1 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ExtSPRClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ExtSS (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ExtSSClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ExtTBR (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Fossilization (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_GeneRate_Dir (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Growth_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_IgrVar (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_IgrBranchRate (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Latent (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Local (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_LocalClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_MixedVar (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_MixedBranchRate (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_MixtureRates (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_MixtureRates_Slider (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_NNI (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_NNIClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_NNI_Hetero (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_NodeSlider (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_NodeSliderClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Nu (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Omega (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Omega_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_OmegaBeta_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_OmegaCat (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_OmegaGamma_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_OmegaM3 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_OmegaNeu (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_OmegaPos (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_OmegaPur (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ParsEraser1 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ParsFossilSPRClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ParsSPR (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ParsSPR1 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ParsSPR2 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ParsSPRClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ParsTBR1 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_ParsTBR2 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_PopSize (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_PosRealLognormal (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_PosRealMultiplier (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_PopSize_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Pinvar (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_RateMult_Dir (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_RateMult_Slider (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_RateShape_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_RealNormal (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_RealSlider (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_RelaxedClockModel (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Revmat_Dir (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Revmat_DirMix (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Revmat_Slider (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Revmat_SplitMerge1 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Revmat_SplitMerge2 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Rho_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Speciation (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Speciation_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Statefreqs (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Statefreqs_Slider (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_StatefreqsSymDirMultistate (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_SwitchRate (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_SwitchRate_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_TK02BranchRate (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_Tratio_Dir (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_TreeStretch (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_TreeLen (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);

#endif  /* __PROPOSAL_H__ */
