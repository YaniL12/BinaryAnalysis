J/A+A/638/A145      GALAH survey. FGK binary stars              (Traven+, 2020)
================================================================================
The GALAH survey: Multiple stars and our Galaxy.
I. A comprehensive method for deriving properties of FGK binary stars.
    Traven G., Feltzing S., Merle T., Van der Swaelmen M., Cotar K., Church R.,
    Zwitter T., Ting Y.-S., Sahlholdt C., Asplund M., Bland-Hawthorn J.,
    De Silva G., Freeman K., Martell S., Sharma S., Zucker D., Buder S.,
    Casey A., D'Orazi V., Kos J., Lewis G., Lin J., Lind K., Simpson J.,
    Stello D., Munari U., Wittenmyer R.A.
    <Astron. Astrophys. 638, A145 (2020)>
    =2020A&A...638A.145T        (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Binaries, spectroscopic ; Photometry, infrared ; Optical ;
              Effective temperatures ; Abundances ; Radial velocities ;
              Stars, diameters
Keywords: methods: data analysis - techniques: radial velocities -
          catalogs - stars: statistics - binaries: spectroscopic

Abstract:
    Binary stellar systems form a large fraction of the Galaxy's stars.
    They are useful as laboratories for studying the physical processes
    taking place within stars, and must be correctly taken into account
    when observations of stars are used to study the structure and
    evolution of the Galaxy. We present a sample of 12760
    well-characterised double-lined spectroscopic binaries that are
    appropriate for statistical studies of the binary populations. They
    were detected as SB2s using a t-distributed stochastic neighbour
    embedding (t-SNE) classification and a cross-correlation analysis of
    GALAH spectra. This sample consists mostly of dwarfs, with a
    significant fraction of evolved stars and several dozen members of the
    giant branch. To compute parameters of the primary and secondary star
    (Teff[1,2], logg[1,2], [Fe/H], Vr[1,2], vmic[1,2], vbroad[1,2],
    R[1,2], and E(B-V)), we used a Bayesian approach that includes a
    parallax prior from Gaia DR2, spectra from GALAH, and apparent
    magnitudes from APASS, Gaia DR2, 2MASS, and WISE. The derived stellar
    properties and their distributions show trends that are expected for a
    population of close binaries (a<10AU) with mass ratios 0.5<=q<=1. The
    derived metallicity of these binary stars is statistically lower than
    that of single dwarf stars from the same magnitude-limited sample.

Description:
    We here analyse a specific data-set: the extended GALAH dataset. This
    consists of stellar spectra from the GALAH survey (reduced as
    explained in Kos et al., 2017MNRAS.464.1259K), apparent magnitudes
    from a variety of photometric catalogues (AAVSO Photometric All Sky
    Survey - APASS; Henden et al. 2016, Cat. II/336, Gaia DR2; Gaia
    Collaboration et al. 2018, Cat. I/345. Two Micron All Sky Survey -
    2MASS; Skrutskie et al. 2006, Cat. VII/233, Wide-field Infrared Survey
    Explorer - WISE; Wright et al. 2010, Cat. II/311), and the parallax
    measurements from Gaia DR2.

    The data provided in this catalogue are described in Table A.1 of the
    paper.

File Summary:
--------------------------------------------------------------------------------
 FileName    Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe          80        .   This file
catalog.dat   2381    12760   Catalogue of 12760 successfully analysed binary
                               stars (described in table A1 of the paper)
--------------------------------------------------------------------------------

See also:
           II/336 : AAVSO Photometric All Sky Survey (APASS) DR9 (Henden+, 2016)
           II/311 : WISE All-Sky Data Release (Cutri+ 2012)
          VII/233 : The 2MASS Extended sources (IPAC/UMass, 2003-2006)
 J/MNRAS/478/4513 : GALAH Survey DR2 (Buder+, 2018)

Byte-by-byte Description of file: catalog.dat
--------------------------------------------------------------------------------
     Bytes Format Units   Label        Explanations
--------------------------------------------------------------------------------
    1-  15  I15   ---     spectID      Unique per-observation star ID
                                        (spectrum ID) (sobject_id)
   17-  32  A16   ---     2MASS        2MASS ID number by default (star_id)
   34-  52  I19   ---     GaiaDR2      Gaia DR2 identifier (gaia_id)
   54-  57  I4    ---     Field        ? GALAH field identification number
                                        (field_id)
   59-  75 F17.13 deg     RAdeg        Right ascension from 2MASS, J2000
                                        (raj2000)
   77-  94 F18.14 deg     DEdeg        Declination from 2MASS, J2000 (dej2000)
   96- 101  F6.3  mag     Bmag         ? B magnitude from APASS (Bmag)
  103- 107  F5.3  mag   e_Bmag         ? Uncertainty of Bmag (e_Bmag)
  109- 114  F6.3  mag     Vmag         ? V magnitude from APASS (Vmag)
  116- 120  F5.3  mag   e_Vmag         ? Uncertainty of Vmag (e_Vmag)
  122- 139 F18.15 mag     GBPmag       ? G_BP magnitude from Gaia DR2 (GBPmag)
  141- 162 F22.20 mag   e_GBPmag       ? Uncertainty of GBPmag (e_GBPmag)
  164- 181 F18.15 mag     Gmag         G magnitude from Gaia DR2 (Gmag)
  183- 204 E22.20 mag   e_Gmag         Uncertainty of Gmag (e_Gmag)
  206- 223 F18.15 mag     GRPmag       ? G_RP magnitude from Gaia DR2 (GRPmag)
  225- 246 F22.20 mag   e_GRPmag       ? Uncertainty of GRPmag (e_GRPmag)
  248- 253  F6.3  mag     Jmag         ? J magnitude from 2MASS (Jmag)
  255- 259  F5.3  mag   e_Jmag         ? Uncertainty of Jmag (e_Jmag)
  261- 266  F6.3  mag     Hmag         ? H magnitude from 2MASS (Hmag)
  268- 272  F5.3  mag   e_Hmag         ? Uncertainty of Hmag (e_Hmag)
  274- 279  F6.3  mag     Kmag         ? Ks magnitude from 2MASS (Kmag)
  281- 285  F5.3  mag   e_Kmag         ? Uncertainty of Kmag (e_Kmag)
  287- 292  F6.3  mag     W1mag        ? W1 magnitude from WISE (W1mag)
  294- 298  F5.3  mag   e_W1mag        ? Uncertainty of W1mag (e_W1mag)
  300- 305  F6.3  mag     W2mag        ? W2 magnitude from WISE (W2mag)
  307- 311  F5.3  mag   e_W2mag        ? Uncertainty of W2mag (e_W2mag)
  313- 333 F21.18 mas     plx           Parallax from Gaia DR2 (parallax)
  335- 354 F20.18 mas   e_plx          Uncertainty of parallax (e_parallax)
  356- 371 F16.12 ---     snc1         ? S/N per pixel in the HERMES
                                        blue channel (sn_c1)
  373- 387 F15.11 ---     snc2         ? S/N per pixel in the HERMES
                                        green channel (sn_c2)
  389- 403 F15.11 ---     snc3         ? S/N per pixel in the HERMES red channel
                                        (sn_c3)
  405- 420 F16.12 ---     snc4         ? S/N per pixel in the HERMES IR channel
                                        (sn_c4)
  422- 423  I2    ---     Flag         [0/32]? Flags in a bitmask format
                                        (flag_cannon)
  425- 442 F18.13 K       Teff1-16     Effective temperature - primary component
                                        16 percentile (teff1_16)
  444- 461 F18.13 K       Teff1-50     Effective temperature - primary component
                                        50 percentile (teff1_50)
  463- 480 F18.13 K       Teff1-84     Effective temperature - primary component
                                        84 percentile (teff1_84)
  482- 499 F18.13 K       Teff1-mean   Effective temperature - primary component
                                       mean value (teff1_mean)
  501- 518 F18.13 K       Teff1-mode   Effective temperature - primary component
                                       mode value (teff1_mode)
  520- 537 F18.13 K       Teff2-16     Effective temperature - secondary
                                       component 16 percentile (teff2_16)
  539- 556 F18.13 K       Teff2-50     Effective temperature - secondary
                                        component 50 percentile (teff2_50)
  558- 575 F18.13 K       Teff2-84     Effective temperature - secondary
                                        component 84 percentile (teff2_84)
  577- 594 F18.13 K       Teff2-mean   Effective temperature - secondary
                                        component mean value (teff2_mean)
  596- 613 F18.13 K       Teff2-mode   Effective temperature - secondary
                                        component mode value (teff2_mode)
  615- 632 F18.16 [cm/s2] logg1-16     Surface gravity - primary component
                                        16 percentile (logg1_16)
  634- 651 F18.16 [cm/s2] logg1-50     Surface gravity - primary component
                                        50 percentile (logg1_50)
  653- 670 F18.16 [cm/s2] logg1-84     Surface gravity - primary component
                                        84 percentile (logg1_84)
  672- 689 F18.16 [cm/s2] logg1-mean   Surface gravity - primary component
                                         mean value (logg1_mean)
  691- 708 F18.16 [cm/s2] logg1-mode   Surface gravity - primary component
                                        mode value (logg1_mode)
  710- 727 F18.16 [cm/s2] logg2-16     Surface gravity - secondary component
                                        16 percentile (logg2_16)
  729- 746 F18.16 [cm/s2] logg2-50     Surface gravity - secondary component
                                         50 percentile (logg2_50)
  748- 765 F18.16 [cm/s2] logg2-84     Surface gravity - secondary component
                                        84 percentile (logg2_84)
  767- 784 F18.16 [cm/s2] logg2-mean   Surface gravity - secondary component
                                         mean value ( logg2_mean)
  786- 803 F18.16 [cm/s2] logg2-mode   Surface gravity - secondary component
                                        mode value (logg2_mode)
  805- 827 E23.20 [-]     [Fe/H]-16    Iron abundance ([Fe/H])
                                        16 percentile (feh_16)
  829- 850 E22.20 [-]     [Fe/H]-50    Iron abundance ([Fe/H])
                                         50 percentile (feh_50)
  852- 874 F23.20 [-]     [Fe/H]-84    Iron abundance ([Fe/H])
                                        84 percentile (feh_84)
  876- 898 F23.20 [-]     [Fe/H]-mean  Iron abundance ([Fe/H])
                                         mean value (feh_mean)
  900- 922 E23.20 [-]     [Fe/H]-mode  Iron abundance ([Fe/H])
                                        mode value (feh_mode)
  924- 948 F25.20 km/s    RV1-16       RV - primary component
                                        16 percentile (V1_16)
  950- 974 F25.20 km/s    RV1-50       RV - primary component
                                         50 percentile (V1_50)
  976- 997 E22.19 km/s    RV1-84       RV - primary component
                                        84 percentile (V1_84)
  999-1020 E22.20 km/s    RV1-mean     RV - primary component
                                        mean value (V1_mean)
 1022-1044 E23.20 km/s    RV1-mode     RV - primary component
                                        mode value (V1_mode)
 1046-1069 F24.19 km/s    RV2-16       RV - secondary component
                                        16 percentile (V2_16)
 1071-1093 F23.18 km/s    RV2-50       RV - secondary component
                                        50 percentile (V2_50)
 1095-1117 F23.18 km/s    RV2-84       RV - secondary component
                                        84 percentile (V2_84)
 1119-1141 F23.18 km/s    RV2-mean     RV - secondary component
                                        mean value (V2_mean)
 1143-1166 F24.19 km/s    RV2-mode     RV - secondary component
                                        mode value (V2_mode)
 1168-1186 F19.17 ---     ratio1-16    Ratio of spectral flux in GALAH band 1
                                        16 percentile (ratio1_16)
 1188-1206 F19.17 ---     ratio1-50    Ratio of spectral flux in GALAH band 1
                                        50 percentile (ratio1_50)
 1208-1226 F19.17 ---     ratio1-84    Ratio of spectral flux in GALAH band 1
                                        84 percentile (ratio1_84)
 1228-1246 F19.17 ---     ratio1-mean  Ratio of spectral flux in GALAH band 1
                                        mean value (ratio1_mean)
 1248-1267 F20.18 ---     ratio1-mode  Ratio of spectral flux in GALAH band 1
                                        mode value (ratio1_mode)
 1269-1288 F20.18 ---     ratio2-16    Ratio of spectral flux in GALAH band 2
                                        16 percentile (ratio2_16)
 1290-1308 F19.17 ---     ratio2-50    Ratio of spectral flux in GALAH band 2
                                        50 percentile (ratio2_50)
 1310-1328 F19.17 ---     ratio2-84    Ratio of spectral flux in GALAH band 2
                                        84 percentile (ratio2_84)
 1330-1348 F19.17 ---     ratio2-mean  Ratio of spectral flux in GALAH band 2
                                        mean value (ratio2_mean)
 1350-1369 F20.18 ---     ratio2-mode  Ratio of spectral flux in GALAH band 2
                                        mode value (ratio2_mode)
 1371-1390 F20.18 ---     ratio3-16    Ratio of spectral flux in GALAH band 3
                                        16 percentile (ratio3_16)
 1392-1410 F19.17 ---     ratio3-50    Ratio of spectral flux in GALAH band 3
                                        50 percentile (ratio3_50)
 1412-1430 F19.17 ---     ratio3-84    Ratio of spectral flux in GALAH band 3
                                        84 percentile (ratio3_84)
 1432-1450 F19.17 ---     ratio3-mean  Ratio of spectral flux in GALAH band 3
                                        mean value (ratio3_mean)
 1452-1470 F19.17 ---     ratio3-mode  Ratio of spectral flux in GALAH band 3
                                        mode value (ratio3_mode)
 1472-1491 F20.18 ---     ratio4-16    Ratio of spectral flux in GALAH band 4
                                        16 percentile (ratio4_16)
 1493-1511 F19.17 ---     ratio4-50    Ratio of spectral flux in GALAH band 4
                                        50 percentile (ratio4_50)
 1513-1531 F19.17 ---     ratio4-84    Ratio of spectral flux in GALAH band 4
                                        84 percentile (ratio4_84)
 1533-1551 F19.17 ---     ratio4-mean  Ratio of spectral flux in GALAH band 4
                                        mean value (ratio4_mean)
 1553-1571 F19.17 ---     ratio4-mode  Ratio of spectral flux in GALAH band 4
                                        mode value (ratio4_mode)
 1573-1590 F18.16 km/s    vmic1-16     Microturbulence - primary component
                                        16 percentile (vmic1_16)
 1592-1609 F18.16 km/s    vmic1-50     Microturbulence - primary component
                                        50 percentile (vmic1_50)
 1611-1628 F18.16 km/s    vmic1-84     Microturbulence - primary component
                                        84 percentile (vmic1_84)
 1630-1647 F18.16 km/s    vmic1-mean   Microturbulence - primary component
                                        mean value (vmic1_mean)
 1649-1666 F18.16 km/s    vmic1-mode   Microturbulence - primary component
                                        mode value (vmic1_mode)
 1668-1685 F18.16 km/s    vmic2-16     Microturbulence - secondary component
                                        16 percentile (vmic2_16)
 1687-1704 F18.16 km/s    vmic2-50     Microturbulence - secondary component
                                        50 percentile (vmic2_50)
 1706-1723 F18.16 km/s    vmic2-84     Microturbulence - secondary component
                                        84 percentile (vmic2_84)
 1725-1742 F18.16 km/s    vmic2-mean   Microturbulence - secondary component
                                        mean value (vmic2_mean)
 1744-1761 F18.16 km/s    vmic2-mode   Microturbulence - secondary component
                                        mode value (vmic2_mode)
 1763-1781 F19.16 km/s    vbroad1-16   Line broadening - primary component
                                        16 percentile (vbroad1_16)
 1783-1801 F19.16 km/s    vbroad1-50   Line broadening - primary component
                                        50 percentile (vbroad1_50)
 1803-1821 F19.16 km/s    vbroad1-84   Line broadening - primary component
                                        84 percentile (vbroad1_84)
 1823-1841 F19.16 km/s    vbroad1-mean Line broadening - primary component
                                        mean value (vbroad1_mean)
 1843-1861 F19.16 km/s    vbroad1-mode Line broadening - primary component
                                        mode value (vbroad1_mode)
 1863-1881 F19.16 km/s    vbroad2-16   Line broadening - secondary component
                                        16 percentile (vbroad2_16)
 1883-1901 F19.16 km/s    vbroad2-50   Line broadening - secondary component
                                        50 percentile (vbroad2_50)
 1903-1921 F19.16 km/s    vbroad2-84   Line broadening - secondary component
                                        84 percentile (vbroad2_84)
 1923-1941 F19.16 km/s    vbroad2-mean Line broadening - secondary component
                                        mean value (vbroad2_mean)
 1943-1961 F19.16 km/s    vbroad2-mode Line broadening - secondary component
                                        mode value (vbroad2_mode)
 1963-1984 E22.20 mag     E(B-V)-16    Interstellar reddening
                                        16 percentile (reddening_16)
 1986-2007 E22.20 mag     E(B-V)-50    Interstellar reddening
                                        50 percentile (reddening_50)
 2009-2030 E22.20 mag     E(B-V)-84    Interstellar reddening
                                        84 percentile (reddening_84)
 2032-2053 E22.20 mag     E(B-V)-mean  Interstellar reddening
                                        mean value (reddening_mean)
 2055-2076 E22.20 mag     E(B-V)-mode  Interstellar reddening
                                        mode value (reddening_mode)
 2078-2098 F21.17 Rsun    R1-16        Stellar radius - primary component
                                        16 percentile (r1_16)
 2100-2119 F20.17 Rsun    R1-50        Stellar radius - primary component
                                        50 percentile (r1_50)
 2121-2141 F21.17 Rsun    R1-84        Stellar radius - primary component
                                        84 percentile (r1_84)
 2143-2163 F21.17 Rsun    R1-mean      Stellar radius - primary component
                                        mean value (r1_mean)
 2165-2185 F21.17 Rsun    R1-mode      Stellar radius - primary component
                                        mode value (r1_mode)
 2187-2207 F21.17 Rsun    R2-16        Stellar radius - secondary component
                                        16 percentile (r2_16)
 2209-2228 F20.17 Rsun    R2-50        Stellar radius - secondary component
                                        50 percentile (r2_50)
 2230-2250 F21.17 Rsun    R2-84        Stellar radius - secondary component
                                        84 percentile (r2_84)
 2252-2272 F21.17 Rsun    R2-mean      Stellar radius - secondary component
                                        mean value (r2_mean)
 2274-2294 F21.17 Rsun    R2-mode      Stellar radius - secondary component
                                        mode value (r2_mode)
 2296-2316 F21.14 ---     chi2         Overall chi2 of the fit
                                        (chi2_binary_pipeline)
 2318-2338 F21.14 ---     chi2sp       Chi2 of the spectroscopic fit
                                        (chi2_binary_pipeline_spec)
 2340-2361 F22.16 ---     chi2ph       Chi2 of the photometric fit
                                        (chi2_binary_pipeline_phot)
 2363-2381 F19.16 ---     RUWE         ? Renormalised unit weight error (ruwe)
--------------------------------------------------------------------------------

Acknowledgements:
    Gregor Traven, gregor.traven(at)astro.lu.se

================================================================================
(End)                                        Patricia Vannier [CDS]  07-Jun-2020
