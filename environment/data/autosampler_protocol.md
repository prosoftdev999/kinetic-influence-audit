# Autosampler carryover qualification

The autosampler qualification immediately preceded the reactor campaign. Its sequence export kept the acquisition positions but lost the vial IDs. `autosampler_nominals.csv` contains the six-channel normalized response expected from each qualification vial when injected after a clean needle. `autosampler_sequence.csv` contains the measured six-channel responses in acquisition order. Every listed vial was injected exactly once.

The qualification model includes one-step needle memory. Channels `polar_1` through `polar_3` share one carryover fraction, and channels `sticky_1` through `sticky_3` share a second fraction. Each fraction is one of the five values 0.18, 0.24, 0.30, 0.36, or 0.42. Except across a deep wash, the predicted response at acquisition i is the current vial's nominal response plus the previous vial's nominal response multiplied by the appropriate carryover fraction. The first acquisition has no previous-vial contribution.

Exactly one programmed deep wash occurred between two adjacent acquisitions. For the acquisition immediately after that wash, the previous-vial contribution is zero. That acquisition then becomes the previous vial for the following acquisition, so ordinary carryover resumes. The independent standard deviation of each normalized channel is 0.003.

Recover the vial order, the two carryover fractions, and the wash location by minimizing

`chi_square = sum(((measured - predicted) / 0.003)^2)`

over all vial permutations, both carryover-grid choices, and every possible wash boundary. Use each vial exactly once. Report the 1-based acquisition number after which the wash occurred. If hypotheses agree within 1e-9 in chi-square, choose the lexicographically smaller `vial_order`; if that is also tied, choose the smaller wash index, then the smaller polar carryover, then the smaller sticky carryover.
