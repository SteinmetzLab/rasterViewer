# Synchronized vs. desynchronized thalamic burst firing: literature notes

Compiled 2026-11-22 for the observation of burst-rich vs. tonic epochs on a deep shank (putative thalamus, no histology yet), with burst-rich epochs apparently splitting into (a) population-synchronized burst volleys at a few Hz and (b) bursts at similar rates but with random relative timing.

**How citations were checked.** Every citation below was checked against its PubMed record (authors, year, title, journal, volume, pages, DOI) via NCBI E-utilities. Unless marked **[full text read]**, the description of each paper comes from its PubMed abstract only. Statements marked **[inference]** are my interpretation, not claims made by the papers.

---

## Short answer

Yes, partly. The literature distinguishes **rhythmic** from **nonrhythmic (arrhythmic)** thalamic bursting, and **synchronized** from **unsynchronized** bursting across neurons. It does not, as far as I found, recognize "synchronized vs. asynchronous burst mode" in awake animals as a named pair of regimes. The closest match to the synchronized regime is the awake-mouse **3-5 Hz "alpha-like" thalamocortical oscillation** (Nestvogel & McCormick 2022). It is carried by rhythmic low-threshold Ca2+ spike (LTS) bursts in dLGN/LP, appears in roughly 1 s episodes after movement stops and as the pupil constricts, and depends on visual cortex. The asynchronous regime fits older behaving-animal reports that spontaneous waking bursts are sparse and "never rhythmic" (Guido & Weyand 1995), and mouse dLGN data showing bursting tracks low or falling arousal on timescales of seconds to minutes (Crombie et al. 2024). Pairwise burst synchrony has been measured directly mainly in anesthetized rat VPm (Whitmire et al. 2016) and in sleeping primates (Liu et al. 2026). For awake mice I did not find a direct measurement; your data could fill that gap.

---

## Most relevant papers

### 1. Nestvogel & McCormick 2022: awake mouse 3-5 Hz alpha-like rhythm driven by rhythmic thalamic LTS bursts
Nestvogel DB, McCormick DA (2022). Visual thalamocortical mechanisms of waking state-dependent activity and alpha oscillations. *Neuron* 110(1):120-138.e4. https://doi.org/10.1016/j.neuron.2021.10.005 (PMID 34687663)

- **Prep:** awake head-fixed mice. Whole-cell recordings in V1 L2/3 and thalamus (6 LP, 2 dLGN cells), plus Neuropixels recordings in dLGN and LP. Mostly no visual stimulation. **[full text read]**
- **What it shows:** The V1 3-5 Hz rhythm is highly coherent with summed thalamic spiking. Each episode opens with a thalamic hyperpolarization followed by 2-7 cycles. Depolarizing phases often trigger LTS bursts, defined as ISI <4 ms after >100 ms silence, which is essentially your criterion. Episodes last about 1 s (median). They cluster after whisking or walking stops, and pupil constriction precedes most of them. Optogenetic silencing of visual cortex abolished the thalamic oscillation, and thalamic muscimol abolished the V1 oscillation. The authors say the role of the thalamic reticular nucleus (TRN) is untested.
- **Relevance:** This is the strongest candidate for your synchronized regime: a population oscillating together at 3-5 Hz, with bursts on each cycle. **[inference]** The paper did not report thalamus-to-thalamus pairwise burst synchrony, and its causal evidence rests on few intracellular thalamic cells.

### 2. Rhythmic vs. nonrhythmic bursting in behaving animals (cat, monkey)
- Guido W, Weyand T (1995). Burst responses in thalamic relay cells of the awake behaving cat. *J Neurophysiol* 74(4):1782-1786. https://doi.org/10.1152/jn.1995.74.4.1782 (PMID 8989413)
- Ramcharan EJ, Gnadt JW, Sherman SM (2000). Burst and tonic firing in thalamic cells of unanesthetized, behaving monkeys. *Vis Neurosci* 17(1):55-62. https://doi.org/10.1017/s0952523800171056 (PMID 10750826)
- Weyand TG, Boudreaux M, Guido W (2001). Burst and tonic response modes in thalamic neurons during sleep and wakefulness. *J Neurophysiol* 85(3):1107-1118. https://doi.org/10.1152/jn.2001.85.3.1107 (PMID 11247981)

- **Prep:** single units in behaving cat LGN (1995, 2001) and behaving macaque LGN and somatosensory thalamus (2000), across wake and sleep.
- **What they show:** In waking cat LGN, spontaneous bursts were infrequent and "never rhythmic" (Guido & Weyand 1995); most waking bursts were visually evoked or followed saccades. Weyand et al. (2001) found bursting during slow-wave sleep came in 1-5 s episodes, about 47% of them rhythmic in the delta range (0.5-4 Hz), and that waking bursts became more likely as vigilance dropped. Ramcharan et al. (2000) found sleep bursting in monkey LGN was *not* rhythmic, while somatosensory thalamus bursting during sleep was "highly rhythmic". In other words, rhythmicity differs between nuclei.
- **Relevance:** These papers established that bursting comes in rhythmic, oscillation-embedded episodes and in isolated, arrhythmic events. That parallels your two regimes, though they analyzed single cells rather than population alignment. **[inference]** The nucleus dependence matters because your shank's nucleus is unknown.

### 3. Whitmire et al. 2016: synchronized bursting across thalamic pairs, gated by adaptation
Whitmire CJ, Waiblinger C, Schwarz C, Stanley GB (2016). Information coding through adaptive gating of synchronized thalamic bursting. *Cell Rep* 14(4):795-807. https://doi.org/10.1016/j.celrep.2015.12.068 (PMID 26776512)

- **Prep:** rat VPm, mostly fentanyl/midazolam/dexmedetomidine anesthesia. Only 7 cells came from awake head-fixed rats. **[full text read]**
- **What it shows:** Sensory adaptation moved VPm along a continuum from burst to tonic firing and lowered pairwise synchrony of burst spikes. Synchrony was a cross-correlogram count within a few ms (12 simultaneously recorded same-barreloid pairs). Burst-rich, synchronized states favored detection; adapted, tonic, less synchronized states favored discrimination. Synchronous bursting appeared mainly with strong stimuli. Burst definition: ISI ~4 ms after 100 ms silence.
- **Relevance:** This gives a direct quantitative framework (pairwise burst-spike synchrony) and a functional reading (detection vs. discrimination). Caveats: the data are stimulus-evoked, anesthetized, and from rat. Spontaneous synchrony was low.

### 4. Liu et al. 2026: burst periodicity and synchrony tested directly across vigilance states
Liu X, Guang J, Israel Z, Wajnsztajn D, Raz A, Bergman H (2026). Entrained cortical delta-spindle activity, not periodic synchrony, prevents arousal by NREM thalamic bursts. *Commun Biol* 9(1):285. https://doi.org/10.1038/s42003-026-09565-3 (PMID 41577949)

- **Prep:** naturally sleeping and waking female macaques (2 animals). Ventral anterior and centromedian thalamus (basal-ganglia-recipient), with polysomnography and field potentials.
- **What it shows (abstract only; full text was not retrievable):** NREM bursting rose, but NREM bursts were "neither periodic nor tightly synchronized". The authors argue that the cortical state at burst onset (delta and spindle context), not burst periodicity or synchrony, explains why NREM bursts do not wake cortex. They frame this as testing a hypothesis that sleep bursts are harmless because they are periodic and synchronized.
- **Relevance:** This is the clearest recent paper treating burst synchrony as something to measure rather than assume. It also warns that "burst-rich" does not imply "synchronized", even in NREM sleep. I could not check their synchrony metric.

### 5. Swadlow & Gusev 2001 and Sherman 2001: the "wake-up call" hypothesis
- Swadlow HA, Gusev AG (2001). The impact of 'bursting' thalamic impulses at a neocortical synapse. *Nat Neurosci* 4(4):402-408. https://doi.org/10.1038/86054 (PMID 11276231)
- Sherman SM (2001). A wake-up call from the thalamus. *Nat Neurosci* 4(4):344-346. https://doi.org/10.1038/85973 (PMID 11276218). This is the News & Views on Swadlow & Gusev, and it is the source of the "wake-up call" phrase.
- Sherman SM (2001). Tonic and burst firing: dual modes of thalamocortical relay. *Trends Neurosci* 24(2):122-126. https://doi.org/10.1016/s0166-2236(00)01714-8 (PMID 11164943)

- **Prep:** awake rabbit VPM to S1 (Swadlow & Gusev); reviews (Sherman).
- **What they show:** In the awake rabbit, a single spontaneous thalamic burst strongly drives its cortical targets. The first spike of a burst is much more effective, and the effect depends on the silent interval before the burst. Sherman's review argues burst mode gives better detectability and tonic mode better linear relay, with corticothalamic and other modulatory inputs switching between them.
- **Relevance:** The wake-up-call idea rests on bursts being *sparse and isolated* (long preceding silence). **[inference]** Your asynchronous regime looks closer to this "detection" use of bursts. Your synchronized regime is the classic low-vigilance oscillatory use, in which cortex is entrained rather than alerted.

### 6. Bezdudnaya et al. 2006 and Bereshpolova et al. 2011: bursting appears within a second of becoming nonalert (awake rabbit LGN)
- Bezdudnaya T, Cano M, Bereshpolova Y, Stoelzel CR, Alonso JM, Swadlow HA (2006). Thalamic burst mode and inattention in the awake LGNd. *Neuron* 49(3):421-432. https://doi.org/10.1016/j.neuron.2006.01.010 (PMID 16446145)
- Bereshpolova Y, Stoelzel CR, Zhuang J, Amitai Y, Alonso JM, Swadlow HA (2011). Getting drowsy? Alert/nonalert transitions and visual thalamocortical network dynamics. *J Neurosci* 31(48):17480-17487. https://doi.org/10.1523/JNEUROSCI.2262-11.2011 (PMID 22131409)

- **Prep:** awake rabbit LGNd (and V1 in 2011), with EEG-defined alert vs. nonalert awake states.
- **What they show:** Within about 1 s of an EEG shift from alert to nonalert, LGN bursting rose sharply and temporal filtering narrowed. In 2011, LGN firing fell and bursting rose in the nonalert state, while L4 spiny-cell rates stayed stable and fast-spiking interneurons became more active.
- **Relevance:** Burst-rich epochs inside wakefulness are linked to EEG-defined inattention, switching on a fast timescale. These papers did not report cross-neuron burst synchrony.

### 7. Crombie et al. 2024: mouse dLGN bursts coupled to pupil dynamics from seconds to minutes
Crombie D, Spacek MA, Leibold C, Busse L (2024). Spiking activity in the visual thalamus is coupled to pupil dynamics across temporal scales. *PLoS Biol* 22(5):e3002614. https://doi.org/10.1371/journal.pbio.3002614 (PMID 38743775)

- **Prep:** awake head-fixed mice on a ball, dLGN, 32-channel silicon probe. Burst definition: >=2 spikes, ISI <4 ms, >=100 ms silence (identical to yours). **[full text read, partially]**
- **What it shows:** Bursts and tonic spikes favor opposite phases of pupil fluctuations at every timescale from about 1 s to about 1000 s. In most pupil-modulated neurons, burst rate peaks at the smallest pupil sizes. Bursts were about 3% of spikes on average. I did not see any analysis of cross-neuron burst synchrony or 3-5 Hz rhythmicity in the portion I read.
- **Relevance:** The closest mouse precedent for your burst-rich vs. tonic epochs lasting tens of seconds to minutes. It also supplies the most obvious covariate (pupil) to test.

### 8. Senzai et al. 2019 and Einstein et al. 2017: the cortical side of the awake mouse 3-6 Hz rhythm
- Senzai Y, Fernandez-Ruiz A, Buzsáki G (2019). Layer-specific physiological features and interlaminar interactions in the primary visual cortex of the mouse. *Neuron* 101(3):500-513.e5. https://doi.org/10.1016/j.neuron.2018.12.009 (PMID 30635232)
- Einstein MC, Polack PO, Tran DT, Golshani P (2017). Visually evoked 3-5 Hz membrane potential oscillations reduce the responsiveness of visual cortex neurons in awake behaving mice. *J Neurosci* 37(20):5084-5098. https://doi.org/10.1523/JNEUROSCI.3868-16.2017 (PMID 28432140)

- **Senzai [full text read]:** frequent 3-6 Hz epochs in V1 during waking immobility, strongest in L4. The rhythm vanishes when the mouse moves and is reduced in NREM sleep. The authors *hypothesize* that it is a precursor of primate alpha. Thalamic origin is assumed, not tested. (The lead said "3-4 Hz"; the paper says 3-6 Hz.)
- **Einstein:** whole-cell recordings in L2/3 of awake mice showed stereotyped 3-5 Hz Vm oscillations lasting 1-2 s, often visually evoked. They reduce responses to preferred stimuli and occur at both high and low arousal (pupil, running).
- **Relevance:** If your synchronized volleys are this rhythm, V1 on another probe should show coincident 3-6 Hz LFP, strongest in L4. **[inference]** Einstein's finding that the rhythm also occurs at high arousal and is visually evoked means pupil alone may not separate the regimes.

### 9. Synchrony of thalamic inputs determines cortical impact
- Alonso JM, Usrey WM, Reid RC (1996). Precisely correlated firing in cells of the lateral geniculate nucleus. *Nature* 383(6603):815-819. https://doi.org/10.1038/383815a0 (PMID 8893005)
- Bruno RM, Sakmann B (2006). Cortex is driven by weak but synchronously active thalamocortical synapses. *Science* 312(5780):1622-1627. https://doi.org/10.1126/science.1124593 (PMID 16778049)
- Wang Q, Webber RM, Stanley GB (2010). Thalamic synchrony and the adaptive gating of information flow to cortex. *Nat Neurosci* 13(12):1534-1541. https://doi.org/10.1038/nn.2670 (PMID 21102447)
- Temereanca S, Brown EN, Simons DJ (2008). Rapid changes in thalamic firing synchrony during repetitive whisker stimulation. *J Neurosci* 28(44):11153-11164. https://doi.org/10.1523/JNEUROSCI.1586-08.2008 (PMID 18971458)

- **Prep:** cat LGN (Alonso; anesthetized per my recollection, not stated in the abstract), rat barrel system in vivo (Bruno & Sakmann; anesthesia details not verified), anesthetized rat VPm and S1 (Wang), rat VPm barreloid pairs (Temereanca).
- **What they show:** Neighboring LGN cells with similar receptive fields often fire within 1 ms of each other, and synchronous spikes drive shared cortical targets more effectively (Alonso). Single thalamocortical synapses are weak, but convergent inputs are numerous and synchronous enough to drive L4 without cortical amplification (Bruno & Sakmann). Adaptation-driven changes in thalamic synchrony track changes in cortical detection vs. discrimination (Wang). Thalamic synchrony varies from millisecond to millisecond with stimulus and adaptation, partly through common input (Temereanca).
- **Relevance:** These papers make the functional stakes concrete. **[inference]** Bursts at the same mean rate should have very different cortical effects in the two regimes, since synchronized burst volleys sum in L4 and asynchronous bursts largely do not. Note that these studies concern millisecond synchrony of single spikes. Burst-onset alignment at the scale of a 3-5 Hz cycle is a coarser kind of synchrony.

### 10. Mechanisms that synchronize thalamic bursting: TRN drive and corticothalamic feedback
- Halassa MM, Siegle JH, Ritt JT, Ting JT, Feng G, Moore CI (2011). Selective optical drive of thalamic reticular nucleus generates thalamic bursts and cortical spindles. *Nat Neurosci* 14(9):1118-1120. https://doi.org/10.1038/nn.2880 (PMID 21785436)
- Contreras D, Destexhe A, Sejnowski TJ, Steriade M (1996). Control of spatiotemporal coherence of a thalamic oscillation by corticothalamic feedback. *Science* 274(5288):771-774. https://doi.org/10.1126/science.274.5288.771 (PMID 8864114)
- Steriade M, McCormick DA, Sejnowski TJ (1993). Thalamocortical oscillations in the sleeping and aroused brain. *Science* 262(5134):679-685. https://doi.org/10.1126/science.8235588 (PMID 8235588)
- Crunelli V, Lőrincz ML, Connelly WM, David F, Hughes SW, Lambert RC, Leresche N, Errington AC (2018). Dual function of thalamic low-vigilance state oscillations: rhythm-regulation and plasticity. *Nat Rev Neurosci* 19(2):107-118. https://doi.org/10.1038/nrn.2017.151 (PMID 29321683)

- **What they show:** Brief optogenetic drive of TRN in behaving mice switched thalamocortical firing from tonic to burst and produced state-dependent spindles (Halassa). In barbiturate-anesthetized cats, removing cortex turned long-range synchronized thalamic spindles into disorganized patterns with low spatiotemporal coherence, while local synchrony between nearby neurons survived (Contreras). Steriade et al. (1993) review how neuromodulatory activation on waking blocks low-frequency synchronized oscillations. Crunelli et al. (2018) review the common thread of thalamic low-vigilance rhythms (<1 Hz to 13 Hz): repetitive Ca2+-spike-driven bursting.
- **Relevance:** **[inference]** Synchronized vs. asynchronous bursting could reflect whether a shared rhythmic inhibitory or excitatory drive is present (TRN, corticothalamic feedback) on top of the hyperpolarized state that permits bursts at all. Neuromodulatory tone would set the hyperpolarization. Contreras et al. imply that synchrony range (local vs. global) is a measurable feature in its own right.

---

## Additional verified papers worth knowing (brief)

- **Higher-order nuclei burst more.** Ramcharan EJ, Gnadt JW, Sherman SM (2005). Higher-order thalamic relays burst more than first-order relays. *PNAS* 102(34):12236-12241. https://doi.org/10.1073/pnas.0502843102 (PMID 16099832). Awake monkey: pulvinar and MD burst far more than LGN, VP, and ventral MGN. Relevant because a deep shank could sit in LP, LD, or other higher-order nuclei.
- **Corticothalamic control of burst vs. tonic mode.** Mease RA, Krieger P, Groh A (2014). *PNAS* 111(18):6798-6803. https://doi.org/10.1073/pnas.1318665111 (L6 stimulation shifts mouse thalamus from bursting to single spikes). Spacek MA et al. (2022). *eLife* 11:e70469. https://doi.org/10.7554/eLife.70469 (CT feedback promotes tonic firing in mouse dLGN; its effects resemble those of behavioral state). Molnár B et al. (2021). *Cereb Cortex Commun* 2(2):tgab020. https://doi.org/10.1093/texcom/tgab020 (LGN TC membrane potential tracks pupil; CT inactivation suppresses the arousal dependence). Reinhold K, Resulaj A, Scanziani M (2023). *J Neurosci* 43(9):1540-1554. https://doi.org/10.1523/JNEUROSCI.2124-21.2022 (state-dependent dLGN responses depend on CT feedback).
- **Arousal and cortical synchrony in quiet wakefulness.** Reimer J et al. (2014). *Neuron* 84(2):355-362. https://doi.org/10.1016/j.neuron.2014.09.033 (pupil constriction goes with more low-frequency oscillation and higher ensemble correlation, switching on a 1-2 s scale). McGinley MJ, David SV, McCormick DA (2015). *Neuron* 87(1):179-192. https://doi.org/10.1016/j.neuron.2015.05.038 (rhythmic slow oscillations at low arousal, stable hyperpolarization at intermediate arousal). Poulet JF, Fernandez LM, Crochet S, Petersen CC (2012). *Nat Neurosci* 15(3):370-372. https://doi.org/10.1038/nn.3035 (cortical desynchronization during active behavior is driven by increased thalamic firing in awake mice).
- **Alpha in cat LGN uses a different burst type.** Hughes SW et al. (2004). *Neuron* 42(2):253-268. https://doi.org/10.1016/s0896-6273(04)00191-6. Lőrincz ML et al. (2009). *Neuron* 63(5):683-696. https://doi.org/10.1016/j.neuron.2009.08.012. Alpha/theta in LGN is driven by *high-threshold* bursts in a subset of TC cells and synchronized by gap junctions. Relay-mode spikes are temporally framed by interneuron inhibition. These bursts are not LTS bursts, and **[inference]** your 100 ms silence criterion would usually exclude them at alpha frequencies.
- **TRN and local sleep.** Fernandez LM et al. (2018). Thalamic reticular control of local sleep in mouse sensory cortex. *eLife* 7:e39111. https://doi.org/10.7554/eLife.39111 (sensory TRN sectors burst more strongly via CaV3.3 and set local NREM spindle vs. delta balance). Vyazovskiy VV et al. (2011). Local sleep in awake rats. *Nature* 472(7344):443-447. https://doi.org/10.1038/nature10009 (cortical neurons go briefly offline with local slow waves after prolonged waking; cortical, not thalamic, recordings).

### Alternative explanations to rule out
- **Head-fixed mice can sleep.** Yüzgeç Ö, Prsa M, Zimmermann R, Huber D (2018). *Curr Biol* 28(3):392-400.e3. https://doi.org/10.1016/j.cub.2017.12.049. Head-fixed mice entered NREM sleep, with pupil size tightly coupled to cortical activity during NREM. **[inference]** Burst-rich epochs lasting tens of seconds to minutes in a 40 min session could be NREM bouts. Synchronized volleys could then be delta or slow-oscillation bursting, and asynchronous bursting could be light NREM, drowsiness, or transitions. This is my top alternative to rule out.
- **Spike-wave discharges.** Letts VA, Beyer BJ, Frankel WN (2014). *Genes Brain Behav* 13(6):519-526. https://doi.org/10.1111/gbb.12142. **[full text read]** 8 of 27 inbred strains had frequent spike-wave discharges (SWDs), typically 6-8 Hz, lasting about 0.6-3 s, mostly during immobility. C57BL/6J was negative. McCafferty C et al. (2018). *Nat Neurosci* 21(5):744-756. https://doi.org/10.1038/s41593-018-0130-4. During absence seizures in behaving animals, TC neurons fired mostly tonically, rarely burst, and were synchronized by cortical drive and TRN feedforward inhibition. **[inference]** SWDs are unlikely if your mice are C57BL/6J, and their frequency (6-8 Hz) is faster than your volleys. Worth confirming strain or background, especially for transgenic crosses.

---

## Testable predictions and analyses for your data

1. **Separate fast synchrony from shared slow modulation.** Compute burst-onset cross-correlograms for all simultaneously recorded pairs, separately within candidate synchronized and asynchronous epochs. Use jitter correction: compare to surrogates that jitter burst onsets by about ±50 ms (keeps slow co-modulation, removes cycle-level alignment) and by about ±500 ms. The synchronized regime should show a central peak above the ±50 ms surrogate plus side peaks at ±200-300 ms. The asynchronous regime should be flat beyond the slow co-modulation.
2. **Population synchrony index on sliding windows.** Bin all burst onsets (for example 25 ms bins), then compute variance/mean against surrogates with independently circularly shifted per-neuron burst trains. Do this in 2-5 s windows. Test whether the index is **bimodal** (two genuine regimes) or **continuous** (one graded variable, as Whitmire et al. found for burst/tonic). A two-state HMM on the index can segment epochs objectively.
3. **Spectrum and autocorrelation of the population burst train.** A 3-5 Hz peak in synchronized epochs supports the Nestvogel/Senzai rhythm. A 0.5-4 Hz peak with long silent periods suggests NREM delta or slow oscillation. A 6-8 Hz peak suggests SWD. Also measure episode durations: Nestvogel reports about 1 s episodes of 2-7 cycles. Sustained volleys over tens of seconds would argue for sleep.
4. **Participation per cycle.** For each population volley, count which neurons burst. Sparse, rotating participation (as in absence seizures and spindles) and all-cell participation imply different mechanisms. Also compare participation against distance along the shank (local vs. global synchrony; Contreras et al. 1996), and across probes if another probe also reaches thalamus.
5. **Check the burst detector's frequency bias.** The 100 ms silence criterion cannot flag more than one burst per cycle above roughly 8-10 Hz, so spindle-frequency bursting would show only the first burst of each spindle. Rerun with relaxed silence (for example 50 ms) and check that the two regimes keep their relative burst rates.
6. **Behavioral and arousal covariates.** Relate the synchronized/asynchronous labels to pupil size and pupil derivative (Crombie: bursts at small or constricting pupil), locomotion and face motion, and time since movement offset (Nestvogel: oscillations follow movement offset and pupil constriction). Prediction: synchronized volleys follow movement offset and pupil constriction, while asynchronous bursting spans broader low-arousal periods. Einstein et al. caution that 3-5 Hz events can also occur at high arousal and can be visually evoked, so check stimulus timing too.
7. **Cortical readout on the other probes.** If you have V1 or other cortex: (a) thalamic volley-triggered cortical LFP and current source density, expecting 3-6 Hz power strongest in L4 for the alpha-like rhythm; (b) cortical DOWN states (population silences of 100-300 ms) and spindles (7-14 Hz), which would indicate NREM; (c) hippocampal sharp-wave ripples or delta if a probe reaches hippocampus. Functional prediction from the synchrony literature: L4 units should respond more per burst during synchronized volleys than during asynchronous bursting at matched burst rates.
8. **Look for TRN-like units on the same or nearby shanks.** TRN cells fire longer, high-rate bursts. **[inference]** If they lead TC volleys in the synchronized regime and are quiet or asynchronous otherwise, that would point to TRN-paced synchrony.
9. **Localize before interpreting.** Higher-order nuclei (LP, LD, POm) burst more than first-order nuclei (Ramcharan 2005) and may differ in rhythmicity (Ramcharan 2000). Histology or atlas registration of the deep shank will change which literature applies.

---

## Leads I could not verify, or that needed correction

- **Senzai et al. 2019, "3-4 Hz".** Corrected: the paper describes a **3-6 Hz** rhythm during waking immobility, largest in L4, which the authors hypothesize is an alpha precursor.
- **Wang, Webber & Stanley 2010.** Verified as Wang **Q**, Webber RM, Stanley GB, *Nat Neurosci* 13:1534 (anesthetized rats).
- **Nestvogel & McCormick, "2022".** Verified as *Neuron* 110(1), January 2022 (online October 2021). It reports rhythmic LTS bursts during the 3-5 Hz rhythm. I did **not** find a direct burst-rate vs. pupil analysis in the text I read; pupil constriction is reported before oscillation onset.
- **Pairwise burst synchrony in LGN pairs from the Alonso, Usrey, or Reid labs.** Not found. Alonso et al. 1996 measures millisecond synchrony of spikes generally, not of bursts specifically. I did not find a paper from these labs quantifying *burst* synchrony between LGN pairs; one may exist that my searches missed.
- **Thalamic bursting vs. pupil from the McCormick, McGinley, or Reimer labs.** Beyond Nestvogel & McCormick 2022 (thalamus) and the cortical papers above (McGinley 2015, Reimer 2014), I found no thalamic burst-vs-pupil paper from these labs. The best mouse thalamic burst-pupil paper I found is from the Busse lab (Crombie et al. 2024).
- **A paper naming an "asynchronous" or "desynchronized burst mode" in awake thalamus.** Not found under that or similar names.
- **Liu et al. 2026 (macaque).** Citation verified. The full text was blocked, so I could not check their synchrony and periodicity metrics or which earlier authors they credit for the "periodic and synchronized" hypothesis.
- **Sherman 2001 (Trends Neurosci).** Verified, but I read only the abstract. I did not confirm whether it explicitly contrasts rhythmic sleep bursting with arrhythmic waking bursting; cite Guido & Weyand 1995, Ramcharan et al. 2000, and Weyand et al. 2001 for that point.
- **Anesthesia details for Bruno & Sakmann 2006 and Alonso et al. 1996.** Not verified from full text; the abstracts do not state them.
- **Whitmire et al. 2016 awake data.** Only 7 awake cells. The text gives the synchrony window as ±5 ms in one place and 10 ms in another; I did not read the supplement.
