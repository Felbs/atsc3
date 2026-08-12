"""ATSC 3.0 L1 signaling syntax, extracted verbatim from the ATSC A/322 standard.

Primary source : ATSC A/322:2026-04 "Physical Layer Protocol", 14 April 2026
                 (lab/spec/A322_2026_tbl.txt, pdftotext -table)
Cross-check    : ATSC A/322:2024-04, 3 April 2024
                 (lab/spec/A322_2024_tbl.txt)

Sections covered
    6.1.2.2   CRC (the CRC-32 used by L1B_crc / L1D_crc)
    6.5.2.1   Table 6.17 -- Ksig for L1-Basic is FIXED at 200 bits
    7.2.3     Table 7.1  -- carrier reduction (Cred_coeff)
    8.5       Table 8.9  -- guard interval durations in samples
    9.2       Table 9.2  -- L1-Basic Signaling Fields and Syntax  (+ Tables 9.3-9.8)
    9.3       Table 9.9  -- L1-Detail Signaling Fields and Syntax (+ Tables 9.10-9.26)

CROSS-EDITION NOTE (2024 -> 2026), the only difference in Table 9.2:
    A/322:2024  ... L1B_first_sub_sbs_last (1), L1B_reserved (48), L1B_crc (32)
    A/322:2026  ... L1B_first_sub_sbs_last (1), L1B_first_sub_mimo_mixed (1),
                    L1B_reserved (47), L1B_crc (32)
    i.e. the 2026 edition carved one bit of the reserved pad out for the new
    L1B_first_sub_mimo_mixed flag.  Every other field name, order and width is
    byte-for-byte identical between the two editions.  Both sum to exactly 200.
    A decoder that treats bit 120 (0-based, MSB-first) as mimo_mixed is correct
    for both editions, because a 2024 transmitter pads reserved with a value the
    receiver is told to ignore -- but note a 2024 transmitter is NOT required to
    set that bit to 0, so only trust it when L1B_version indicates the 2026
    structure.

    L1-Detail changed more between editions:
      * L1D_version: 1 (2024) -> 2 (2026)
      * 2026 appends a SECOND per-subframe loop after L1D_bsid carrying
        L1D_mimo_mixed and the per-PLP L1D_plp_mimo block.  A 2024-era parser
        sees that whole loop as part of L1D_reserved, which is exactly the
        forward-compatibility mechanism described in 9.3.1.

All bit widths below are transcribed from the "No. of Bits" column.
Page headers/footers in the source text ("ATSC A/322:2026-04",
"Physical Layer Protocol", "14 April 2026", bare page numbers) were stripped.
"""

# ---------------------------------------------------------------------------
# Section 9.2 -- Table 9.2  L1-Basic Signaling Fields and Syntax
# ---------------------------------------------------------------------------
#
# The one conditional in L1-Basic is on L1B_frame_length_mode.  Both arms are
# 23 bits wide, so the total is a constant 200 bits either way:
#
#     if (L1B_frame_length_mode == 0) {          # time-aligned frames
#         L1B_frame_length               10
#         L1B_excess_samples_per_symbol  13
#     } else {                                   # symbol-aligned frames
#         L1B_time_offset                16
#         L1B_additional_samples          7
#     }

#: Structured form: entries are ("field", name, bits) or
#: ("if_frame_length_mode", {0: [...], 1: [...]}).
L1_BASIC_STRUCT = [
    ("field", "L1B_version", 3),
    ("field", "L1B_mimo_scattered_pilot_encoding", 1),
    ("field", "L1B_lls_flag", 1),
    ("field", "L1B_time_info_flag", 2),
    ("field", "L1B_return_channel_flag", 1),
    ("field", "L1B_papr_reduction", 2),
    ("field", "L1B_frame_length_mode", 1),
    ("if_frame_length_mode", {
        0: [("field", "L1B_frame_length", 10),
            ("field", "L1B_excess_samples_per_symbol", 13)],
        1: [("field", "L1B_time_offset", 16),
            ("field", "L1B_additional_samples", 7)],
    }),
    ("field", "L1B_num_subframes", 8),
    ("field", "L1B_preamble_num_symbols", 3),
    ("field", "L1B_preamble_reduced_carriers", 3),
    ("field", "L1B_L1_Detail_content_tag", 2),
    ("field", "L1B_L1_Detail_size_bytes", 13),
    ("field", "L1B_L1_Detail_fec_type", 3),
    ("field", "L1B_L1_Detail_additional_parity_mode", 2),
    ("field", "L1B_L1_Detail_total_cells", 19),
    ("field", "L1B_first_sub_mimo", 1),
    ("field", "L1B_first_sub_miso", 2),
    ("field", "L1B_first_sub_fft_size", 2),
    ("field", "L1B_first_sub_reduced_carriers", 3),
    ("field", "L1B_first_sub_guard_interval", 4),
    ("field", "L1B_first_sub_num_ofdm_symbols", 11),
    ("field", "L1B_first_sub_scattered_pilot_pattern", 5),
    ("field", "L1B_first_sub_scattered_pilot_boost", 3),
    ("field", "L1B_first_sub_sbs_first", 1),
    ("field", "L1B_first_sub_sbs_last", 1),
    ("field", "L1B_first_sub_mimo_mixed", 1),   # 2026 only; part of L1B_reserved in 2024
    ("field", "L1B_reserved", 47),              # 48 in A/322:2024-04
    ("field", "L1B_crc", 32),
]


def _flatten_basic(frame_length_mode):
    out = []
    for ent in L1_BASIC_STRUCT:
        if ent[0] == "field":
            out.append((ent[1], ent[2]))
        else:
            for sub in ent[1][frame_length_mode]:
                out.append((sub[1], sub[2]))
    return out


#: Flat transmission order for the time-aligned case (L1B_frame_length_mode = 0).
#: This is the usual real-world US broadcast configuration.
L1_BASIC_FIELDS = _flatten_basic(0)

#: Flat transmission order for the symbol-aligned case (L1B_frame_length_mode = 1).
L1_BASIC_FIELDS_SYMBOL_ALIGNED = _flatten_basic(1)

#: A/322:2024-04 variant (no L1B_first_sub_mimo_mixed; L1B_reserved is 48 bits).
L1_BASIC_FIELDS_2024 = [
    (name, 48 if name == "L1B_reserved" else bits)
    for (name, bits) in L1_BASIC_FIELDS
    if name != "L1B_first_sub_mimo_mixed"
]

#: Section 6.5.2.1 / Table 6.17: "The value of Ksig for L1-Basic is fixed as 200".
#: Section 9.1 preamble text: "The length of L1-Basic signaling is fixed at 200 bits."
L1_BASIC_KSIG_BITS = 200


# ---------------------------------------------------------------------------
# L1-Basic semantics (Tables 9.3-9.8, plus the shared Tables 9.11-9.16)
# ---------------------------------------------------------------------------

#: Table 8.9 -- guard interval durations in samples, keyed by GI pattern name.
GI_PATTERN_SAMPLES = {
    "GI1_192": 192,
    "GI2_384": 384,
    "GI3_512": 512,
    "GI4_768": 768,
    "GI5_1024": 1024,
    "GI6_1536": 1536,
    "GI7_2048": 2048,
    "GI8_2432": 2432,    # not available for 8K FFT
    "GI9_3072": 3072,    # not available for 8K FFT
    "GI10_3648": 3648,   # not available for 8K FFT
    "GI11_4096": 4096,   # not available for 8K FFT
    "GI12_4864": 4864,   # 32K FFT only
}

#: FFT sizes for which each GI pattern is legal (check marks in Table 8.9).
GI_PATTERN_FFT_ALLOWED = {
    "GI1_192": ("8K", "16K", "32K"),
    "GI2_384": ("8K", "16K", "32K"),
    "GI3_512": ("8K", "16K", "32K"),
    "GI4_768": ("8K", "16K", "32K"),
    "GI5_1024": ("8K", "16K", "32K"),
    "GI6_1536": ("8K", "16K", "32K"),
    "GI7_2048": ("8K", "16K", "32K"),
    "GI8_2432": ("16K", "32K"),
    "GI9_3072": ("16K", "32K"),
    "GI10_3648": ("16K", "32K"),
    "GI11_4096": ("16K", "32K"),
    "GI12_4864": ("32K",),
}

#: Table 9.13 -- index -> GI pattern name (shared by L1B_first_sub_guard_interval
#: and L1D_guard_interval).
GUARD_INTERVAL_MAP = {
    0b0000: "Reserved",
    0b0001: "GI1_192",
    0b0010: "GI2_384",
    0b0011: "GI3_512",
    0b0100: "GI4_768",
    0b0101: "GI5_1024",
    0b0110: "GI6_1536",
    0b0111: "GI7_2048",
    0b1000: "GI8_2432",
    0b1001: "GI9_3072",
    0b1010: "GI10_3648",
    0b1011: "GI11_4096",
    0b1100: "GI12_4864",
    0b1101: "Reserved",
    0b1110: "Reserved",
    0b1111: "Reserved",
}

#: Convenience: index -> guard interval length in samples (None where reserved).
GUARD_INTERVAL_SAMPLES = {
    k: GI_PATTERN_SAMPLES.get(v) for k, v in GUARD_INTERVAL_MAP.items()
}

#: Table 9.14 -- scattered pilot pattern for SISO (shared by
#: L1B_first_sub_scattered_pilot_pattern and L1D_scattered_pilot_pattern).
SCATTERED_PILOT_PATTERN_SISO = {
    0b00000: "SP3_2",
    0b00001: "SP3_4",
    0b00010: "SP4_2",
    0b00011: "SP4_4",
    0b00100: "SP6_2",
    0b00101: "SP6_4",
    0b00110: "SP8_2",
    0b00111: "SP8_4",
    0b01000: "SP12_2",
    0b01001: "SP12_4",
    0b01010: "SP16_2",
    0b01011: "SP16_4",
    0b01100: "SP24_2",
    0b01101: "SP24_4",
    0b01110: "SP32_2",
    0b01111: "SP32_4",
}
SCATTERED_PILOT_PATTERN_SISO.update({v: "Reserved" for v in range(0b10000, 0b100000)})

#: Table 9.15 -- scattered pilot pattern for MIMO (same indices, "MP" names).
SCATTERED_PILOT_PATTERN_MIMO = {
    k: v.replace("SP", "MP") for k, v in SCATTERED_PILOT_PATTERN_SISO.items()
    if v != "Reserved"
}
SCATTERED_PILOT_PATTERN_MIMO.update({v: "Reserved" for v in range(0b10000, 0b100000)})

#: Derived from the SPx_y naming: (DX, DY) of the scattered pilot lattice.
SCATTERED_PILOT_DX_DY = {
    k: (int(v[2:].split("_")[0]), int(v.split("_")[1]))
    for k, v in SCATTERED_PILOT_PATTERN_SISO.items() if v != "Reserved"
}

#: Section 7.2.3 / Table 7.1 -- reduced_carriers index is Cred_coeff (0..4);
#: 5,6,7 are reserved.  NoC = NoCmax - Cred_coeff * Cunit,
#: Cunit = 96 (8K), 192 (16K), 384 (32K).
CARRIER_REDUCTION_UNIT = {"8K": 96, "16K": 192, "32K": 384}
NOC_MAX = {"8K": 6913, "16K": 13825, "32K": 27649}   # Cred_coeff = 0
NUM_CARRIERS = {          # Table 7.1: Cred_coeff -> NoC per FFT size
    0: {"8K": 6913, "16K": 13825, "32K": 27649},
    1: {"8K": 6817, "16K": 13633, "32K": 27265},
    2: {"8K": 6721, "16K": 13441, "32K": 26881},
    3: {"8K": 6625, "16K": 13249, "32K": 26497},
    4: {"8K": 6529, "16K": 13057, "32K": 26113},
}
REDUCED_CARRIERS_MAP = {
    0: "Cred_coeff=0 (no reduction)",
    1: "Cred_coeff=1",
    2: "Cred_coeff=2",
    3: "Cred_coeff=3",
    4: "Cred_coeff=4",
    5: "Reserved",
    6: "Reserved",
    7: "Reserved",
}

#: Table 9.12 -- FFT size (shared by L1B_first_sub_fft_size and L1D_fft_size).
FFT_SIZE_MAP = {0b00: "8K", 0b01: "16K", 0b10: "32K", 0b11: "Reserved"}
FFT_SIZE_POINTS = {0b00: 8192, 0b01: 16384, 0b10: 32768, 0b11: None}

#: Table 9.11 -- MISO option (shared by L1B_first_sub_miso and L1D_miso).
MISO_MAP = {
    0b00: "No MISO",
    0b01: "MISO with 64 coefficients",
    0b10: "MISO with 256 coefficients",
    0b11: "Reserved",
}

#: Table 9.16 -- scattered pilot boost, power in dB, indexed
#: [SP pattern name][boost value].  Values 101/110/111 are RFU.
SCATTERED_PILOT_BOOST_DB = {
    "SP3_2":  {0: 0.00, 1: 0.00, 2: 1.40, 3: 2.20, 4: 2.90},
    "SP3_4":  {0: 0.00, 1: 1.40, 2: 2.90, 3: 3.80, 4: 4.40},
    "SP4_2":  {0: 0.00, 1: 0.60, 2: 2.10, 3: 3.00, 4: 3.60},
    "SP4_4":  {0: 0.00, 1: 2.10, 2: 3.60, 3: 4.40, 4: 5.10},
    "SP6_2":  {0: 0.00, 1: 1.60, 2: 3.10, 3: 4.00, 4: 4.60},
    "SP6_4":  {0: 0.00, 1: 3.00, 2: 4.50, 3: 5.40, 4: 6.00},
    "SP8_2":  {0: 0.00, 1: 2.20, 2: 3.80, 3: 4.60, 4: 5.30},
    "SP8_4":  {0: 0.00, 1: 3.60, 2: 5.10, 3: 6.00, 4: 6.60},
    "SP12_2": {0: 0.00, 1: 3.20, 2: 4.70, 3: 5.60, 4: 6.20},
    "SP12_4": {0: 0.00, 1: 4.50, 2: 6.00, 3: 6.90, 4: 7.50},
    "SP16_2": {0: 0.00, 1: 3.80, 2: 5.30, 3: 6.20, 4: 6.80},
    "SP16_4": {0: 0.00, 1: 5.20, 2: 6.70, 3: 7.60, 4: 8.20},
    "SP24_2": {0: 0.00, 1: 4.70, 2: 6.20, 3: 7.10, 4: 7.70},
    "SP24_4": {0: 0.00, 1: 6.10, 2: 7.60, 3: 8.50, 4: 9.10},
    "SP32_2": {0: 0.00, 1: 5.40, 2: 6.90, 3: 7.70, 4: 8.40},
    "SP32_4": {0: 0.00, 1: 6.70, 2: 8.20, 3: 9.10, 4: 9.70},
}
# The MIMO patterns MPx_y share the same boost values as SPx_y (Table 9.16 lists
# them as "SP3_2 / MP3_2" etc.).
SCATTERED_PILOT_BOOST_DB.update(
    {k.replace("SP", "MP"): v for k, v in list(SCATTERED_PILOT_BOOST_DB.items())}
)

L1_BASIC_SEMANTICS = {
    # Table 9.3
    "L1B_mimo_scattered_pilot_encoding": {
        0: "Walsh-Hadamard pilots or no MIMO Subframes",
        1: "Null pilots",
    },
    "L1B_lls_flag": {
        0: "No LLS signaling in the current Frame",
        1: "LLS signaling carried in this Frame",
    },
    # Table 9.4
    "L1B_time_info_flag": {
        0b00: "Time information is not included in the current Frame",
        0b01: "Time information included, signaled to ms precision",
        0b10: "Time information included, signaled to us precision",
        0b11: "Time information included, signaled to ns precision",
    },
    "L1B_return_channel_flag": {
        0: "Dedicated return channel (DRC) not supported",
        1: "Dedicated return channel (DRC) supported",
    },
    # Table 9.5
    "L1B_papr_reduction": {
        0b00: "No PAPR reduction used",
        0b01: "Tone reservation only",
        0b10: "ACE only",
        0b11: "Both TR and ACE",
    },
    "L1B_frame_length_mode": {
        0: "Time-aligned Frame (excess samples distributed to guard intervals); "
           "L1B_frame_length + L1B_excess_samples_per_symbol present",
        1: "Symbol-aligned Frame (no excess sample distribution); "
           "L1B_time_offset + L1B_additional_samples present",
    },
    # Table 9.6
    "L1B_L1_Detail_fec_type": {
        0b000: "Mode 1",
        0b001: "Mode 2",
        0b010: "Mode 3",
        0b011: "Mode 4",
        0b100: "Mode 5",
        0b101: "Mode 6",
        0b110: "Mode 7",
        0b111: "Reserved",
    },
    # Table 9.7
    "L1B_L1_Detail_additional_parity_mode": {
        0b00: "K = 0 (no additional parity used)",
        0b01: "K = 1",
        0b10: "K = 2",
        0b11: "Reserved for future use",
    },
    # Table 9.8 (in conjunction with L1B_first_sub_mimo_mixed)
    "L1B_first_sub_mimo": {
        0: "First subframe includes one or more PLPs without MIMO processing",
        1: "MIMO processing applied to all PLPs in the first subframe",
    },
    "L1B_first_sub_mimo_mixed": {
        0: "All PLPs in the subframe either use MIMO or all do not",
        1: "PLPs using and not using MIMO are multiplexed within the subframe "
           "(requires LDM; MIMO on Enhanced PLPs only)",
    },
    "L1B_first_sub_miso": MISO_MAP,
    "L1B_first_sub_fft_size": FFT_SIZE_MAP,
    "L1B_preamble_reduced_carriers": REDUCED_CARRIERS_MAP,
    "L1B_first_sub_reduced_carriers": REDUCED_CARRIERS_MAP,
    "L1B_first_sub_guard_interval": GUARD_INTERVAL_MAP,
    "L1B_first_sub_scattered_pilot_pattern": SCATTERED_PILOT_PATTERN_SISO,
    "L1B_first_sub_scattered_pilot_boost": {
        0b000: "boost index 0 (see Table 9.16, pattern-dependent)",
        0b001: "boost index 1",
        0b010: "boost index 2",
        0b011: "boost index 3",
        0b100: "boost index 4",
        0b101: "Reserved for future use",
        0b110: "Reserved for future use",
        0b111: "Reserved for future use",
    },
    "L1B_first_sub_sbs_first": {
        0: "First symbol of the first Subframe is not a Subframe boundary symbol",
        1: "First symbol of the first Subframe is a Subframe boundary symbol",
    },
    "L1B_first_sub_sbs_last": {
        0: "Last symbol of the first Subframe is not a Subframe boundary symbol",
        1: "Last symbol of the first Subframe is a Subframe boundary symbol",
    },
}

#: Fields whose signaled value is "one less than" the real quantity (9.2.x text).
L1_BASIC_OFF_BY_ONE = {
    "L1B_num_subframes": "actual number of Subframes = value + 1",
    "L1B_preamble_num_symbols": "actual number of Preamble OFDM symbols = value + 1",
    "L1B_first_sub_num_ofdm_symbols":
        "actual number of data payload OFDM symbols (incl. SBS) = value + 1",
}

#: Other scaling rules quoted from 9.2.1 / 9.2.2.
L1_BASIC_UNITS = {
    "L1B_frame_length": "units of 5 ms; frame length = value * 5 ms; 10 <= value <= 1000",
    "L1B_excess_samples_per_symbol": "extra samples added to the GI of each "
                                     "non-Preamble OFDM symbol",
    "L1B_time_offset": "sample periods (at the frame BSR) from the preceding/coincident "
                       "millisecond boundary to the leading edge of the Frame",
    "L1B_additional_samples": "additional samples at end of Frame; shall be 0 in this "
                              "version of the specification",
    "L1B_L1_Detail_size_bytes": "size of L1-Detail information in bytes, excluding "
                                "additional parity; minimum 25",
    "L1B_L1_Detail_total_cells": "total OFDM cells of coded+modulated L1-Detail for this "
                                 "Frame plus modulated additional parity for the next",
    "L1B_version": "shall be 1 for this version of the specification",
}


# ---------------------------------------------------------------------------
# Section 9.3 -- Table 9.9  L1-Detail Signaling Fields and Syntax
# ---------------------------------------------------------------------------
#
# Node forms used below:
#   ("field", name, bits)
#   ("if",   condition_str, [then...])
#   ("if",   condition_str, [then...], [else...])
#   ("for",  loop_var, count_expr_str, [body...])
#
# Loop counts follow the spec's inclusive "for (i = 0 .. N)" convention, i.e.
# the body runs N+1 times.  Condition/count strings are transcribed from the
# spec and are meant to be evaluated against already-parsed field values.

L1_DETAIL_FIELDS = [
    ("field", "L1D_version", 4),
    ("field", "L1D_num_rf", 3),
    ("for", "L1D_rf_id", "1 .. L1D_num_rf", [          # runs L1D_num_rf times
        ("field", "L1D_bonded_bsid", 16),
        ("field", "reserved", 3),
    ]),
    ("if", "L1B_time_info_flag != 0b00", [
        ("field", "L1D_time_sec", 32),
        ("field", "L1D_time_msec", 10),
        ("if", "L1B_time_info_flag != 0b01", [
            ("field", "L1D_time_usec", 10),
            ("if", "L1B_time_info_flag != 0b10", [
                ("field", "L1D_time_nsec", 10),
            ]),
        ]),
    ]),
    ("for", "i", "0 .. L1B_num_subframes", [           # runs L1B_num_subframes+1 times
        ("if", "i > 0", [
            ("field", "L1D_mimo", 1),
            ("field", "L1D_miso", 2),
            ("field", "L1D_fft_size", 2),
            ("field", "L1D_reduced_carriers", 3),
            ("field", "L1D_guard_interval", 4),
            ("field", "L1D_num_ofdm_symbols", 11),
            ("field", "L1D_scattered_pilot_pattern", 5),
            ("field", "L1D_scattered_pilot_boost", 3),
            ("field", "L1D_sbs_first", 1),
            ("field", "L1D_sbs_last", 1),
        ]),
        ("if", "L1B_num_subframes > 0", [
            ("field", "L1D_subframe_multiplex", 1),
        ]),
        ("field", "L1D_frequency_interleaver", 1),
        ("if", "(i == 0 and (L1B_first_sub_sbs_first or L1B_first_sub_sbs_last)) or "
               "(i > 0 and (L1D_sbs_first or L1D_sbs_last))", [
            ("field", "L1D_sbs_null_cells", 13),
        ]),
        ("field", "L1D_num_plp", 6),
        ("for", "j", "0 .. L1D_num_plp", [             # runs L1D_num_plp+1 times
            ("field", "L1D_plp_id", 6),
            ("field", "L1D_plp_lls_flag", 1),
            ("field", "L1D_plp_layer", 2),
            ("field", "L1D_plp_start", 24),
            ("field", "L1D_plp_size", 24),
            ("field", "L1D_plp_scrambler_type", 2),
            ("field", "L1D_plp_fec_type", 4),
            ("if", "L1D_plp_fec_type in (0,1,2,3,4,5)", [
                ("field", "L1D_plp_mod", 4),
                ("field", "L1D_plp_cod", 4),
            ]),
            ("field", "L1D_plp_TI_mode", 2),
            ("if", "L1D_plp_TI_mode == 0b00", [
                ("field", "L1D_plp_fec_block_start", 15),
            ], [
                ("if", "L1D_plp_TI_mode == 0b01", [
                    ("field", "L1D_plp_CTI_fec_block_start", 22),
                ]),
            ]),
            ("if", "L1D_num_rf > 0", [
                ("field", "L1D_plp_num_channel_bonded", 3),
                ("if", "L1D_plp_num_channel_bonded > 0", [
                    ("field", "L1D_plp_channel_bonding_format", 2),
                    ("for", "k", "0 .. L1D_plp_num_channel_bonded", [
                        ("field", "L1D_plp_bonded_rf_id", 3),
                    ]),
                ]),
            ]),
            ("if", "(i == 0 and L1B_first_sub_mimo == 1) or (i > 0 and L1D_mimo == 1)", [
                ("field", "L1D_plp_mimo_stream_combining", 1),
                ("field", "L1D_plp_mimo_IQ_interleaving", 1),
                ("field", "L1D_plp_mimo_PH", 1),
            ]),
            ("if", "L1D_plp_layer == 0", [
                ("field", "L1D_plp_type", 1),
                ("if", "L1D_plp_type == 1", [
                    ("field", "L1D_plp_num_subslices", 14),
                    ("field", "L1D_plp_subslice_interval", 24),
                ]),
                ("if", "L1D_plp_TI_mode in (0b01, 0b10) and L1D_plp_mod == 0b0000", [
                    ("field", "L1D_plp_TI_extended_interleaving", 1),
                ]),
                ("if", "L1D_plp_TI_mode == 0b01", [
                    ("field", "L1D_plp_CTI_depth", 3),
                    ("field", "L1D_plp_CTI_start_row", 11),
                ], [
                    ("if", "L1D_plp_TI_mode == 0b10", [
                        ("field", "L1D_plp_HTI_inter_subframe", 1),
                        ("field", "L1D_plp_HTI_num_ti_blocks", 4),
                        ("field", "L1D_plp_HTI_num_fec_blocks_max", 12),
                        ("if", "L1D_plp_HTI_inter_subframe == 0", [
                            ("field", "L1D_plp_HTI_num_fec_blocks", 12),
                        ], [
                            ("for", "k", "0 .. L1D_plp_HTI_num_ti_blocks", [
                                ("field", "L1D_plp_HTI_num_fec_blocks", 12),
                            ]),
                        ]),
                        ("field", "L1D_plp_HTI_cell_interleaver", 1),
                    ]),
                ]),
            ], [
                ("field", "L1D_plp_ldm_injection_level", 5),
            ]),
        ]),
    ]),
    ("field", "L1D_bsid", 16),
    # ---- The block below exists only in A/322:2026-04 (L1D_version = 2). ----
    ("for", "i", "0 .. L1B_num_subframes", [
        ("if", "i > 0", [
            ("field", "L1D_mimo_mixed", 1),
        ]),
        ("if", "(i == 0 and L1B_first_sub_mimo_mixed == 1) or "
               "(i > 0 and L1D_mimo_mixed == 1)", [
            ("for", "j", "0 .. L1D_num_plp", [
                ("field", "L1D_plp_mimo", 1),
                ("if", "L1D_plp_mimo == 1", [
                    ("field", "L1D_plp_mimo_stream_combining", 1),
                    ("field", "L1D_plp_mimo_IQ_interleaving", 1),
                    ("field", "L1D_plp_mimo_PH", 1),
                ]),
            ]),
        ]),
    ]),
    # ------------------------------------------------------------------------
    ("field", "L1D_reserved", "as needed"),   # pads to 8*L1B_L1_Detail_size_bytes
    ("field", "L1D_crc", 32),
]

#: Same structure as A/322:2024-04 (L1D_version = 1): identical up to L1D_bsid,
#: then straight to L1D_reserved / L1D_crc with no mimo_mixed loop.
L1_DETAIL_FIELDS_2024 = [
    n for n in L1_DETAIL_FIELDS
    if not (n[0] == "for" and n[1] == "i" and any(
        s[0] == "if" and "mimo_mixed" in s[1] for s in n[3]))
]


L1_DETAIL_SEMANTICS = {
    "L1D_version": {1: "A/322:2024-04 structure", 2: "A/322:2026-04 structure"},
    "L1D_num_rf": {
        0: "Channel bonding not used for the current Frame",
        1: "Bonded with one other channel (max for this version of the spec)",
    },
    # Table 9.10 (with L1D_mimo)
    "L1D_mimo": {
        0: "Subframe includes one or more PLPs without MIMO processing",
        1: "MIMO processing applied to all PLPs in the Subframe",
    },
    "L1D_mimo_mixed": {
        0: "All PLPs in the Subframe either use MIMO or all do not",
        1: "PLPs using and not using MIMO are multiplexed within the Subframe",
    },
    # Table 9.11 / 9.12 / 9.13 / 9.14 -- shared with the L1B first_sub_* fields
    "L1D_miso": MISO_MAP,
    "L1D_fft_size": FFT_SIZE_MAP,
    "L1D_reduced_carriers": REDUCED_CARRIERS_MAP,
    "L1D_guard_interval": GUARD_INTERVAL_MAP,
    "L1D_scattered_pilot_pattern": SCATTERED_PILOT_PATTERN_SISO,
    "L1D_scattered_pilot_boost": {
        0b000: "boost index 0 (see Table 9.16, pattern-dependent)",
        0b001: "boost index 1",
        0b010: "boost index 2",
        0b011: "boost index 3",
        0b100: "boost index 4",
        0b101: "Reserved for future use",
        0b110: "Reserved for future use",
        0b111: "Reserved for future use",
    },
    "L1D_sbs_first": {0: "First symbol is not a Subframe boundary symbol",
                      1: "First symbol is a Subframe boundary symbol"},
    "L1D_sbs_last": {0: "Last symbol is not a Subframe boundary symbol",
                     1: "Last symbol is a Subframe boundary symbol"},
    "L1D_subframe_multiplex": {
        0: "Subframe is TDM / concatenated in time with adjacent Subframes",
        1: "Reserved for future use",
    },
    "L1D_frequency_interleaver": {0: "Frequency interleaver bypassed",
                                  1: "Frequency interleaver enabled"},
    # Table 9.18
    "L1D_plp_scrambler_type": {
        0b00: "Scrambler defined in Section 5.2.3",
        0b01: "Reserved for future use",
        0b10: "Reserved for future use",
        0b11: "Reserved for future use",
    },
    # Table 9.19
    "L1D_plp_fec_type": {
        0b0000: "BCH + 16K LDPC",
        0b0001: "BCH + 64K LDPC",
        0b0010: "CRC + 16K LDPC",
        0b0011: "CRC + 64K LDPC",
        0b0100: "16K LDPC only",
        0b0101: "64K LDPC only",
        **{v: "Reserved for future use" for v in range(0b0110, 0b10000)},
    },
    # Table 9.20 (SISO)
    "L1D_plp_mod": {
        0b0000: "QPSK",
        0b0001: "16QAM-NUC",
        0b0010: "64QAM-NUC",
        0b0011: "256QAM-NUC",
        0b0100: "1024QAM-NUC",
        0b0101: "4096QAM-NUC",
        **{v: "Reserved" for v in range(0b0110, 0b10000)},
    },
    # Table 9.22
    "L1D_plp_cod": {
        0b0000: "2/15",
        0b0001: "3/15",
        0b0010: "4/15",
        0b0011: "5/15",
        0b0100: "6/15",
        0b0101: "7/15",
        0b0110: "8/15",
        0b0111: "9/15",
        0b1000: "10/15",
        0b1001: "11/15",
        0b1010: "12/15",
        0b1011: "13/15",
        **{v: "Reserved" for v in range(0b1100, 0b10000)},
    },
    # Table 9.23
    "L1D_plp_TI_mode": {
        0b00: "No time interleaving (neither CTI nor HTI)",
        0b01: "Convolutional time interleaving (CTI) mode",
        0b10: "Hybrid time interleaving (HTI) mode",
        0b11: "Reserved for future use",
    },
    # Table 9.25
    "L1D_plp_channel_bonding_format": {
        0b00: "Plain channel bonding",
        0b01: "SNR averaged channel bonding",
        0b10: "Reserved for future use",
        0b11: "Reserved for future use",
    },
    # Table 9.26 -- Nrows for the convolutional time interleaver
    "L1D_plp_CTI_depth": {
        0b000: "512",
        0b001: "724",
        0b010: "887 (non-extended interleaving) or 1254 (extended interleaving)",
        0b011: "1024 (non-extended interleaving) or 1448 (extended interleaving)",
        0b100: "Reserved for future use",
        0b101: "Reserved for future use",
        0b110: "Reserved for future use",
        0b111: "Reserved for future use",
    },
    "L1D_plp_layer": {
        0: "Core Layer",
        1: "Enhanced Layer (LDM applied to this Subframe)",
        2: "Reserved in this version of the specification",
        3: "Reserved in this version of the specification",
    },
    "L1D_plp_type": {
        0: "Non-dispersed PLP (contiguous cells, no subslicing)",
        1: "Dispersed PLP (subslicing used)",
    },
    "L1D_plp_lls_flag": {
        0: "PLP does not carry LLS information",
        1: "PLP does carry LLS information (not necessarily in this Frame)",
    },
    "L1D_plp_mimo": {0: "MIMO not used for this PLP", 1: "MIMO used for this PLP"},
    "L1D_plp_mimo_stream_combining": {0: "Stream combining not used",
                                      1: "Stream combining used"},
    "L1D_plp_mimo_IQ_interleaving": {0: "IQ polarization interleaving not used",
                                     1: "IQ polarization interleaving used"},
    "L1D_plp_mimo_PH": {0: "Phase hopping not used", 1: "Phase hopping used"},
    "L1D_plp_TI_extended_interleaving": {0: "Extended interleaving not used",
                                         1: "Extended interleaving used"},
    "L1D_plp_HTI_inter_subframe": {0: "Intra-subframe HTI", 1: "Inter-subframe HTI"},
    # Table 9.24 -- LDM injection level in dB
    "L1D_plp_ldm_injection_level": {
        0b00000: 0.0,  0b00001: 0.5,  0b00010: 1.0,  0b00011: 1.5,
        0b00100: 2.0,  0b00101: 2.5,  0b00110: 3.0,  0b00111: 3.5,
        0b01000: 4.0,  0b01001: 4.5,  0b01010: 5.0,  0b01011: 6.0,
        0b01100: 7.0,  0b01101: 8.0,  0b01110: 9.0,  0b01111: 10.0,
        0b10000: 11.0, 0b10001: 12.0, 0b10010: 13.0, 0b10011: 14.0,
        0b10100: 15.0, 0b10101: 16.0, 0b10110: 17.0, 0b10111: 18.0,
        0b11000: 19.0, 0b11001: 20.0, 0b11010: 21.0, 0b11011: 22.0,
        0b11100: 23.0, 0b11101: 24.0, 0b11110: 25.0, 0b11111: "Reserved",
    },
}

#: Off-by-one / unit rules from Section 9.3.x.
L1_DETAIL_OFF_BY_ONE = {
    "L1D_num_plp": "actual number of PLPs in the Subframe = value + 1",
    "L1D_num_ofdm_symbols": "actual number of data payload OFDM symbols (incl. SBS) "
                            "= value + 1",
    "L1D_plp_num_subslices": "actual number of subslices = value + 1 "
                             "(value 0 is reserved)",
    "L1D_plp_HTI_num_ti_blocks": "number of TI blocks = value + 1",
    "L1D_plp_HTI_num_fec_blocks": "number of FEC blocks = value + 1",
    "L1D_plp_HTI_num_fec_blocks_max": "max number of FEC blocks = value + 1",
}
# NOTE: there is no L1D_num_subframes field.  The subframe count comes from
# L1B_num_subframes in L1-Basic (subframe count = L1B_num_subframes + 1).


# ---------------------------------------------------------------------------
# Section 6.1.2.2 -- CRC (quoted)
# ---------------------------------------------------------------------------
CRC32_SPEC_QUOTE = """\
A/322, Section 6.1.2.2 "CRC":

  "When a CRC is used for the Outer Code, a 32-bit CRC shall be used. The CRC
   shall be computed as illustrated in Figure 6.4 and shall implement a feedback
   shift register characterized by the CRC code polynomial. The generator
   polynomial of degree n Gcrc(x) can be expressed as:
       Gcrc(x) = x^n + g_{n-1} x^{n-1} + g_{n-2} x^{n-2} + ..... + g2 x^2 + g1 x + 1"

  "At the beginning of the computation (before the first data bit is input) all
   register stage contents shall be initialized to one. After applying the first
   bit (MSB first) of the data block to the input, the shift clock causes the
   register to shift its contents by one stage towards b_{n-1} while loading the
   tapped stages with the result of the appropriate operations. After the last
   data bit of the block is input, the contents of the register stages are read
   out to provide the 32 CRC bits {b_i, i = 0,1...31) that shall then be appended
   to the data prior to inner encoding. The appended bits shall be ordered from
   the most significant bit (b31) to the least significant bit (b0). For the
   CRC-32 used, all values of gi = 0 except for: g21, g16, g11 which shall have a
   value of one. Thus the actual generator polynomial shall be:
       Gcrc(x) = x^32 + x^21 + x^16 + x^11 + 1"

L1B_crc (Section 9.2.4):
  "This field shall contain the CRC value as computed according to Section
   6.1.2.2 over the contents of L1-Basic excluding the L1B_crc field."

L1D_crc (Section 9.3.1):
  "This field shall contain the CRC value as computed in Section 6.1.2.2 over the
   contents of L1-Detail excluding the L1D_crc field."
"""

#: Gcrc(x) = x^32 + x^21 + x^16 + x^11 + 1, without the implicit x^32 term.
#: (This is NOT the common CRC-32/MPEG-2 polynomial 0x04C11DB7 -- A/322 uses its
#: own sparse polynomial.  Init = all ones; MSB-first; no final XOR stated.)
CRC32_POLY = (1 << 32) | (1 << 21) | (1 << 16) | (1 << 11) | 1
CRC32_POLY_TRUNCATED = 0x00210801   # bits 21, 16, 11, 0
CRC32_WIDTH = 32
CRC32_INIT = 0xFFFFFFFF
CRC32_REFLECT_IN = False
CRC32_REFLECT_OUT = False
CRC32_XOR_OUT = 0x00000000


def atsc3_l1_crc32(bits):
    """Compute the A/322 6.1.2.2 CRC-32 over an MSB-first iterable of 0/1 ints.

    Direct transcription of the shift-register description: register initialised
    to all ones, data applied MSB first, no final inversion.
    """
    reg = CRC32_INIT
    for b in bits:
        fb = ((reg >> 31) & 1) ^ (b & 1)
        reg = ((reg << 1) & 0xFFFFFFFF)
        if fb:
            reg ^= CRC32_POLY_TRUNCATED
    return reg


# ---------------------------------------------------------------------------
# Section 6.5.2.1 -- Table 6.17, L1-Basic FEC modes (useful for demod)
# ---------------------------------------------------------------------------
L1_BASIC_FEC_MODES = {
    # mode -> (Ksig, code length, code rate, constellation, length in cells)
    1: (200, 16200, "3/15 (Type A)", "QPSK", 3820),
    2: (200, 16200, "3/15 (Type A)", "QPSK", 934),
    3: (200, 16200, "3/15 (Type A)", "QPSK", 484),
    4: (200, 16200, "3/15 (Type A)", "NUC_16_8/15", 259),
    5: (200, 16200, "3/15 (Type A)", "NUC_64_9/15", 163),
    6: (200, 16200, "3/15 (Type A)", "NUC_256_9/15", 112),
    7: (200, 16200, "3/15 (Type A)", "NUC_256_13/15", 69),
}


# ---------------------------------------------------------------------------
# Self-check
# ---------------------------------------------------------------------------

def l1_basic_total_bits(fields=None):
    """Total width of the L1-Basic field list.  Must equal 200."""
    return sum(w for _, w in (fields if fields is not None else L1_BASIC_FIELDS))


def _check():
    ok = True

    total0 = l1_basic_total_bits(L1_BASIC_FIELDS)
    total1 = l1_basic_total_bits(L1_BASIC_FIELDS_SYMBOL_ALIGNED)
    total24 = l1_basic_total_bits(L1_BASIC_FIELDS_2024)

    print("A/322 L1 signaling syntax -- self check")
    print("=" * 66)
    print("L1-Basic fields (frame_length_mode=0, time-aligned): %d fields"
          % len(L1_BASIC_FIELDS))
    print("L1-Basic total bits (mode 0, 2026 edition) = %d  (Ksig = %d)  %s"
          % (total0, L1_BASIC_KSIG_BITS,
             "OK" if total0 == L1_BASIC_KSIG_BITS else "MISMATCH"))
    print("L1-Basic total bits (mode 1, 2026 edition) = %d  %s"
          % (total1, "OK" if total1 == L1_BASIC_KSIG_BITS else "MISMATCH"))
    print("L1-Basic total bits (A/322:2024-04)        = %d  %s"
          % (total24, "OK" if total24 == L1_BASIC_KSIG_BITS else "MISMATCH"))
    ok &= (total0 == total1 == total24 == L1_BASIC_KSIG_BITS)

    # widths must be small positive ints
    bad = [(n, w) for n, w in L1_BASIC_FIELDS if not (isinstance(w, int) and 1 <= w <= 64)]
    print("Bit-width sanity (1..64, int): %s" % ("OK" if not bad else "BAD %r" % bad))
    ok &= not bad

    # guard interval menu
    published_gi = {192, 384, 512, 768, 1024, 1536, 2048, 2432, 3072, 3648, 4096, 4864}
    got_gi = {v for v in GUARD_INTERVAL_SAMPLES.values() if v is not None}
    print("Guard-interval menu: %s" % ("OK" if got_gi == published_gi
                                       else "BAD %r" % (got_gi ^ published_gi)))
    ok &= got_gi == published_gi

    # scattered pilot menu
    published_sp = {"SP3_2", "SP3_4", "SP4_2", "SP4_4", "SP6_2", "SP6_4",
                    "SP8_2", "SP8_4", "SP12_2", "SP12_4", "SP16_2", "SP16_4",
                    "SP24_2", "SP24_4", "SP32_2", "SP32_4"}
    got_sp = {v for v in SCATTERED_PILOT_PATTERN_SISO.values() if v != "Reserved"}
    print("Scattered-pilot menu: %s" % ("OK" if got_sp == published_sp
                                        else "BAD %r" % (got_sp ^ published_sp)))
    ok &= got_sp == published_sp

    # semantics keys must be real field names
    names = {n for n, _ in L1_BASIC_FIELDS} | {"L1B_time_offset",
                                               "L1B_additional_samples"}
    unknown = sorted(set(L1_BASIC_SEMANTICS) - names)
    print("Semantics keys resolve to fields: %s"
          % ("OK" if not unknown else "UNKNOWN %r" % unknown))
    ok &= not unknown

    # Table 7.1 consistency: NoC = NoCmax - Cred_coeff * Cunit
    derived_ok = all(
        NUM_CARRIERS[c][f] == NOC_MAX[f] - c * CARRIER_REDUCTION_UNIT[f]
        for c in range(5) for f in ("8K", "16K", "32K"))
    print("Table 7.1 vs NoC = NoCmax - Cred*Cunit: %s" % ("OK" if derived_ok else "BAD"))
    ok &= derived_ok

    print("=" * 66)
    print("RESULT: %s" % ("ALL CHECKS PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    print()
    for _name, _bits in L1_BASIC_FIELDS:
        print("  %-42s %2d" % (_name, _bits))
    print()
    sys.exit(_check())
