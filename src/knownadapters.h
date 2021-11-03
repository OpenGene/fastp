#ifndef KNOWN_ADAPTERS_H
#define KNOWN_ADAPTERS_H

// some adapter sequences are from
// https://github.com/stephenturner/adapters/blob/master/adapters_combined_256_unique.fasta

#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string>

inline map<string, string> getKnownAdapter() {
  map<string, string> knownAdapters;

  knownAdapters["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"] =
      ">Illumina TruSeq Adapter Read 1";
  knownAdapters["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"] =
      ">Illumina TruSeq Adapter Read 2";
  knownAdapters["GATCGTCGGACTGTAGAACTCTGAACGTGTAGA"] =
      ">Illumina Small RNA Adapter Read 2";
  knownAdapters["AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"] =
      ">Illumina DpnII expression PCR Primer 2 | >Illumina NlaIII expression "
      "PCR Primer 2 | >Illumina Small RNA PCR Primer 2 | >Illumina DpnII Gex "
      "PCR Primer 2 | >Illumina NlaIII Gex PCR Primer 2";
  knownAdapters["AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"] =
      ">Illumina RNA PCR Primer";
  knownAdapters["AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"] =
      ">TruSeq_Universal_Adapter | >PrefixPE/1 | >PCR_Primer1 | >Illumina "
      "Single End PCR Primer 1 | >Illumina Paried End PCR Primer 1 | >Illumina "
      "Multiplexing PCR Primer 1.01 | >TruSeq Universal Adapter | "
      ">TruSeq_Universal_Adapter | >PrefixPE/1 | >PCR_Primer1";
  knownAdapters["AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTAGAT"
                "CGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"] =
      ">pcr_dimer";
  knownAdapters["AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCAAG"
                "CAGAAGACGGCATACGAGCTCTTCCGATCT"] = ">PCR_Primers";
  knownAdapters["ACACTCTTTCCCTACACGACGCTCTTCCGATCT"] =
      ">Illumina Single End Sequencing Primer | >Illumina Paired End Adapter 1 "
      "| >Illumina Paried End Sequencing Primer 1 | >Illumina Multiplexing "
      "Adapter 2 | >Illumina Multiplexing Read1 Sequencing Primer";
  knownAdapters["AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"] =
      ">PE2_rc | >TruSeq3_IndexedAdapter | >PE2_rc | >TruSeq3_IndexedAdapter";
  knownAdapters
      ["AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">Reverse_adapter";
  knownAdapters["AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"] = ">TruSeq2_PE_r";
  knownAdapters
      ["AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">PCR_Primer2_rc";
  knownAdapters
      ["AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAA"] =
          ">PhiX_read1_adapter";
  knownAdapters["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"] =
      ">PE1_rc | >TruSeq3_UniversalAdapter | >PE1_rc | "
      ">TruSeq3_UniversalAdapter";
  knownAdapters["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">PCR_Primer1_rc";
  knownAdapters
      ["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA"] =
          ">PhiX_read2_adapter";
  knownAdapters["AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"] = ">TruSeq2_SE";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 35";
  knownAdapters["CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC"] =
      ">Illumina PCR Primer Index 10";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 10";
  knownAdapters["CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC"] =
      ">Illumina PCR Primer Index 2";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 2";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 38";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 27";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 25";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 31";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 44";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 37";
  knownAdapters["CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC"] =
      ">Illumina PCR Primer Index 6";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 6";
  knownAdapters["CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC"] =
      ">Illumina PCR Primer Index 5";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 5";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 23";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 30";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 21";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 42";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 33";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"] =
          ">PrefixPE/2 | >PCR_Primer2 | >Illumina Paired End PCR Primer 2 | "
          ">PrefixPE/2 | >PCR_Primer2";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 22";
  knownAdapters["CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC"] =
      ">Illumina PCR Primer Index 1";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 1";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 17";
  knownAdapters["CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC"] =
      ">Illumina PCR Primer Index 9";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 9";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 47";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 28";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 45";
  knownAdapters["CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC"] =
      ">Illumina PCR Primer Index 7";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 7";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 34";
  knownAdapters["CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC"] =
      ">Illumina PCR Primer Index 3";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 3";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 18";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 24";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 26";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 43";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 14";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 16";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 20";
  knownAdapters["CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC"] =
      ">Illumina PCR Primer Index 11";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 11";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 39";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 41";
  knownAdapters["CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC"] =
      ">Illumina PCR Primer Index 12";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 12";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 29";
  knownAdapters["CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC"] =
      ">Illumina PCR Primer Index 8";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 8";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 46";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 40";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 15";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 32";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 48";
  knownAdapters["CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC"] =
      ">Illumina PCR Primer Index 4";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 4";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 36";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 13";
  knownAdapters
      ["CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"] =
          ">RNA PCR Primer, Index 19";
  knownAdapters["CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"] =
      ">Illumina Single End Adapter 2 | >Illumina Single End PCR Primer 2";
  knownAdapters["CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT"] =
      ">ABI Solid3 Adapter B";
  knownAdapters["CCGACAGGTTCAGAGTTCTACAGTCCGACATG"] =
      ">Illumina NlaIII expression Sequencing Primer | >Illumina NlaIII Gex "
      "Sequencing Primer";
  knownAdapters["CCGAGCCCACGAGACAAGAGGCAATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N711 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N711 | "
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N711 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N711";
  knownAdapters["CCGAGCCCACGAGACACTCGCTAATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N716";
  knownAdapters["CCGAGCCCACGAGACACTGAGCGATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N724";
  knownAdapters["CCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N703 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N703 | "
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N703 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N703";
  knownAdapters["CCGAGCCCACGAGACATCTCAGGATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N715";
  knownAdapters["CCGAGCCCACGAGACATGCGCAGATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N722";
  knownAdapters["CCGAGCCCACGAGACCAGAGAGGATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N708";
  knownAdapters["CCGAGCCCACGAGACCCTAAGACATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N726";
  knownAdapters["CCGAGCCCACGAGACCGAGGCTGATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N710 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N710 | "
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N710 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N710";
  knownAdapters["CCGAGCCCACGAGACCGATCAGTATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N727";
  knownAdapters["CCGAGCCCACGAGACCGGAGCCTATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N720";
  knownAdapters["CCGAGCCCACGAGACCGTACTAGATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N702 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N702 | "
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N702 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N702";
  knownAdapters["CCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N707 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N707 | "
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N707 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N707";
  knownAdapters["CCGAGCCCACGAGACGCGTAGTAATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N719";
  knownAdapters["CCGAGCCCACGAGACGCTACGCTATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N709";
  knownAdapters["CCGAGCCCACGAGACGCTCATGAATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N714";
  knownAdapters["CCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N705 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N705 | "
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N705 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N705";
  knownAdapters["CCGAGCCCACGAGACGGAGCTACATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N718";
  knownAdapters["CCGAGCCCACGAGACGTAGAGGAATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N712 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N712 | "
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N712 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N712";
  knownAdapters["CCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N701 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N701 | "
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N701 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N701";
  knownAdapters["CCGAGCCCACGAGACTACGCTGCATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N721";
  knownAdapters["CCGAGCCCACGAGACTAGCGCTCATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N723";
  knownAdapters["CCGAGCCCACGAGACTAGGCATGATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N706 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N706 | "
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N706 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N706";
  knownAdapters["CCGAGCCCACGAGACTCCTGAGCATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N704 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N704 | "
      ">I7_Primer_Nextera_XT_and_Nextera_Enrichment_N704 | "
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N704";
  knownAdapters["CCGAGCCCACGAGACTCGACGTCATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N729";
  knownAdapters["CCGAGCCCACGAGACTGCAGCTAATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Primer_Nextera_XT_Index_Kit_v2_N728";
  knownAdapters["CGACAGGTTCAGAGTTCTACAGTCCGACGATC"] =
      ">Illumina DpnII expression Sequencing Primer | >Illumina Small RNA "
      "Sequencing Primer | >Illumina DpnII Gex Sequencing Primer";
  knownAdapters["CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"] =
      ">Illumina Paired End Sequencing Primer 2";
  knownAdapters["CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"] =
      ">Clontech Universal Primer Mix Long";
  knownAdapters["CTGAGCGGGCTGGCAAGGCAGACCGATCTCGTATGCCGTCTTCTGCTTG"] =
      ">I7_Adapter_Nextera_No_Barcode";
  knownAdapters["CTGATGGCGCGAGGGAGGCGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Adapter_Nextera";
  knownAdapters["CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG"] = ">ABI Solid3 Adapter A";
  knownAdapters["CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"] =
      ">I7_Nextera_Transposase_1 | >Trans2_rc | >I7_Nextera_Transposase_1 | "
      ">Trans2_rc";
  knownAdapters["CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGC"] =
      ">I7_Nextera_Transposase_2";
  knownAdapters["CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC"] =
      ">I5_Nextera_Transposase_2";
  knownAdapters["CTGTCTCTTATACACATCTGACGCTGCCGACGA"] =
      ">I5_Nextera_Transposase_1 | >Trans1_rc | >I5_Nextera_Transposase_1 | "
      ">Trans1_rc";
  knownAdapters["GACGCTGCCGACGAACTCTAGGGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S516";
  knownAdapters["GACGCTGCCGACGAAGAGGATAGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]503 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S503 | "
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]503 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S503";
  knownAdapters["GACGCTGCCGACGAAGCTAGAAGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S515";
  knownAdapters["GACGCTGCCGACGAAGGCTTAGGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]508 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S508 | "
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]508 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S508";
  knownAdapters["GACGCTGCCGACGAATAGAGAGGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]502 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S502 | "
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]502 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S502";
  knownAdapters["GACGCTGCCGACGAATAGCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S520";
  knownAdapters["GACGCTGCCGACGAATTAGACGGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S510";
  knownAdapters["GACGCTGCCGACGACGGAGAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S511";
  knownAdapters["GACGCTGCCGACGACTAGTCGAGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S513";
  knownAdapters["GACGCTGCCGACGACTCCTTACGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]505 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S505 | "
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]505 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S505";
  knownAdapters["GACGCTGCCGACGACTTAATAGGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S518";
  knownAdapters["GACGCTGCCGACGAGCGATCTAGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]501";
  knownAdapters["GACGCTGCCGACGATAAGGCTCGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S521";
  knownAdapters["GACGCTGCCGACGATACTCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]507 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S507 | "
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]507 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S507";
  knownAdapters["GACGCTGCCGACGATATGCAGTGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]506 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S506 | "
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]506 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S506";
  knownAdapters["GACGCTGCCGACGATCGCATAAGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S522";
  knownAdapters["GACGCTGCCGACGATCTACTCTGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]504";
  knownAdapters["GACGCTGCCGACGATCTTACGCGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]517 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S517 | "
      ">I5_Primer_Nextera_XT_and_Nextera_Enrichment_[N/S/E]517 | "
      ">I5_Primer_Nextera_XT_Index_Kit_v2_S517";
  knownAdapters["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"] =
      ">Nextera_LMP_Read1_External_Adapter | >Illumina Multiplexing Index "
      "Sequencing Primer";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_5 | >TruSeq Adapter, Index 5";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_25";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq Adapter, Index 25";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_8 | >TruSeq Adapter, Index 8";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_13";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq Adapter, Index 13";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_14";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq Adapter, Index 14";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_1_6 | >TruSeq Adapter, Index 1";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_15";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq Adapter, Index 15";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_27";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq Adapter, Index 27";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_7 | >TruSeq Adapter, Index 7";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCACTCTTCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq Adapter, Index 23";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_16";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq Adapter, Index 16";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_2 | >TruSeq Adapter, Index 2";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_22";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq Adapter, Index 22";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_12 | >TruSeq Adapter, Index 12";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_23";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_9 | >TruSeq Adapter, Index 9";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_6 | >TruSeq Adapter, Index 6";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_11 | >TruSeq Adapter, Index 11";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_18_7";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq Adapter, Index 18";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_19";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq Adapter, Index 19";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_20";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq Adapter, Index 20";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_21";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq Adapter, Index 21";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_10 | >TruSeq Adapter, Index 10";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_4 | >TruSeq Adapter, Index 4";
  knownAdapters
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">TruSeq_Adapter_Index_3 | >TruSeq Adapter, Index 3";
  knownAdapters["GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"] =
      ">Illumina Paired End Adapter 2";
  knownAdapters["GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"] =
      ">Nextera_LMP_Read2_External_Adapter";
  knownAdapters["GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"] =
      ">Illumina Single End Adapter 1";
  knownAdapters["GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"] = ">Trans2";
  knownAdapters["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"] =
      ">PrefixPE/2 | >PE2 | >Illumina Multiplexing PCR Primer 2.01 | >Illumina "
      "Multiplexing Read2 Sequencing Primer | >PrefixPE/2 | >PE2";
  knownAdapters["TACACTCTTTCCCTACACGACGCTCTTCCGATCT"] =
      ">PrefixPE/1 | >PE1 | >PrefixPE/1 | >PE1";
  knownAdapters["TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"] =
      ">RNA_PCR_Primer_(RP1)_part_#_15013198";
  knownAdapters["TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"] = ">Trans1";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_5_(RPI5)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_25_(RPI25)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_8_(RPI8)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTCAAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_13_(RPI13)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTTCCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_14_(RPI14)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_1_(RPI1)_2,9";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGAGCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_26_(RPI26)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGTCAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_15_(RPI15)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATTCCTATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_27_(RPI27)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAAAAGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_28_(RPI28)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAACTAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_29_(RPI29)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACCGGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_30_(RPI30)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACGATATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_31_(RPI31)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACTCAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_32_(RPI32)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_7_(RPI7)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGGCGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_33_(RPI33)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATGGCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_34_(RPI34)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATTTTATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_35_(RPI35)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCAACAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_36_(RPI36)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCGTCCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_16_(RPI16)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_2_(RPI2)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGGAATATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_37_(RPI37)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGTACGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_22_(RPI22)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTAGCTATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_38_(RPI38)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTATACATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_39_(RPI39)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTCAGAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_40_(RPI40)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_12_(RPI12)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_41_(RPI41)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGAGTGGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_23_(RPI23)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_9_(RPI9)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_6_(RPI6)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_11_(RPI11)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGTAGCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_24_(RPI24)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTAGAGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_17_(RPI17)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_18_(RPI18)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGAAAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_19_(RPI19)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGGCCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_20_(RPI20)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTTTCGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_21_(RPI21)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAATCGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_42_(RPI42)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTACAGCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_43_(RPI43)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_10_(RPI10)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTATAATATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_44_(RPI44)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCATTCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_45_(RPI45)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCCCGAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_46_(RPI46)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGAAGATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_47_(RPI47)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGGCAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_48_(RPI48)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_4_(RPI4)";
  knownAdapters
      ["TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"] =
          ">RNA_PCR_Primer_Index_3_(RPI3)";
  knownAdapters["TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC"] = ">FlowCell1";
  knownAdapters["TTTTTTTTTTCAAGCAGAAGACGGCATACGA"] = ">FlowCell2";
  knownAdapters["AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"] =
      ">MGI/BGI adapter (forward)";
  knownAdapters["AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"] =
      ">MGI/BGI adapter (reverse)";

  return knownAdapters;
}
#endif