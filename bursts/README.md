# Burst detection: choosing a definition (2026-11-22)

Goal: classify every spike as burst or tonic, to show them separately (or overlaid in two
colors) in the viewer. This folder compares three definitions on FD_016. The viewer now has
these modes (its "Spikes" menu), using the LTS definition by default with adjustable
parameters; its detection matches `detect_bursts.py` exactly (checked by the self-test).

- `detect_bursts.py`: the detector, and its self-tests on synthetic spike trains.
- `burst_diagnostics.py`: makes `figs/` and `summary.json` (about 50 s).

## Definitions compared

A burst is a run of at least K spikes of one unit, every ISI at most M, whose first spike
follows at least S of silence from that unit. All spikes in the run are burst spikes; the first
is the burst onset. Everything else is tonic.

| Name | Silence S | Max ISI M | Min spikes K |
|---|---|---|---|
| LTS, the standard thalamic convention | >= 100 ms | <= 4 ms | 2 |
| Relaxed LTS | >= 50 ms | <= 6 ms | 2 |
| No silence (the first proposal) | none | <= 8 ms | 3 |

The LTS convention is from Lu, Guido & Sherman 1992 (J Neurophysiol 68:2185). That paper says
a short ISI alone was not enough, and that low-threshold bursts were also marked by a prior
silent period. Later studies state the criterion as ISIs under 4 ms after more than 100 ms of
silence: Reinagel & Reid 2000 (cat LGN), and Whitmire et al. 2016 (Cell Reports, awake rat VPm).
A review by Alitto & Usrey gives 50 ms of silence. I could not read the 1992 methods section to
confirm its exact silence value.

The silence requirement is physiological, not cosmetic. T-type Ca2+ channels must be held
hyperpolarized for tens of ms to de-inactivate, so a low-threshold burst comes after a pause.

## What the figures show

1. **ISI return maps** (`figs/1_isi_return_maps.png`).
   - Probe 6 shank 0 has a separate cluster of spikes with a long ISI before (about 50 ms to
     10 s) and a short ISI after (about 2 to 4 ms). These are burst onsets: the classic LTS
     signature.
   - Units more than 2000 µm above the tip have no such cluster; their ISIs form one continuous
     diagonal.
   - The cluster's following ISIs reach about 5 ms, so a 4 ms cut trims its upper edge
     slightly.
2. **Burst fraction vs position** (`figs/2_burst_fraction_vs_depth.png`).
   - Under the LTS definition, Probe 6 shank 0 has many units between 450 and 1700 µm with 20
     to 85% of their spikes in bursts.
   - Strongly bursting units also appear near the tips of Probes 3, 4 and 5 (below 400 µm),
     probably also thalamic.
   - **Not only thalamus:** some superficial units on Probes 1 and 2 (about 2300 to 2600 µm)
     also reach 30 to 80%. These may be intrinsically bursting cortical cells (layer 5). The
     definition describes a spike pattern; it cannot tell a T-type burst from a cortical burst.
3. **Example raster** (`figs/3_example_raster_P6S0.png`), the window of your screenshot,
   1110.1 to 1111.3 s.
   - All three definitions pick out the burst episodes from about 1110.65 s, at 1200 to
     1750 µm.
   - The no-silence definition adds a few events in the 500 to 900 µm band before the
     transition. It also misses the 2-spike bursts.
4. **Burst structure** (`figs/4_burst_structure.png`).
   - On Probe 6 shank 0, LTS bursts are 42% doublets and 58% three or more spikes. Their ISIs
     are about 2.3 ms and lengthen slightly within the burst (2.33, 2.37, 2.63 ms medians), as
     LTS bursts do.
   - Superficial units' "bursts" are 80% doublets.
   - Under the no-silence definition, superficial units' runs have ISIs of about 4.8 ms that do
     not lengthen: fast tonic firing, not bursts.
5. **Time course** (`figs/5_timecourse_P6S0_deep.png`), Probe 6 shank 0 below 1400 µm.
   - Over the whole session, firing alternates between burst-rich and tonic epochs, lasting
     tens of seconds to minutes.
   - 1050 to 1095 s is burst-rich, then 1100 to 1160 s is mostly tonic. The 1110.6 s burst
     episode in your screenshot is a brief one inside a tonic epoch.
6. **Sensitivity** (`figs/6_sensitivity.png`): % of spikes in bursts on a grid of silence x max
   ISI.
   - Bursting is most specific to deep units (deep / superficial ratio 2.0 to 2.9) with
     ISI <= 3 or 4 ms. With 6 or 8 ms it drops to 1.4 to 2.0.
   - The silence requirement matters less to that ratio than the max ISI does. But without it,
     deep units' burst fraction triples (4.3% to 13.5% at 4 ms), mostly from fast tonic runs.
   - Depth is only a rough stand-in for brain area here.

## Recommendation

Start with the LTS convention (100 ms, 4 ms, 2 spikes) and make the three numbers adjustable in
the viewer. Your first proposal (no silence, 8 ms, at least 3 spikes) misses doublets, which are
42% of LTS bursts on Probe 6 shank 0. It also counts fast tonic runs as bursts.

## Caveats

- No histology: brain areas are inferred from position along the shank and firing patterns.
- Spike sorting: the later spikes of a burst are often smaller, so a sorter can split a burst
  across units, or drop its later spikes. That would undercount spikes per burst and miss some
  bursts.
- The first spike of each unit's recording is treated as preceded by silence. Negligible.
