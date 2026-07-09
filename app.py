"""Streamlit app for NGS_LibraryQC."""

from __future__ import annotations

import gzip
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from ngs_libraryqc import clean_seq, run_counting


st.set_page_config(page_title="NGS_LibraryQC", layout="wide")
st.title("NGS_LibraryQC")
st.caption("Anchor-based NGS library QC for pooled amplicon libraries")


def save_uploaded_file(uploaded_file) -> str:
    suffix = ".fastq.gz" if uploaded_file.name.endswith(".gz") else ".fastq"
    temp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    temp.write(uploaded_file.getvalue())
    temp.close()
    return temp.name


def save_uploaded_csv(uploaded_file) -> str:
    temp = tempfile.NamedTemporaryFile(delete=False, suffix=".csv")
    temp.write(uploaded_file.getvalue())
    temp.close()
    return temp.name


with st.sidebar:
    st.header("Input FASTQ")
    fastq_files = st.file_uploader(
        "Upload FASTQ or FASTQ.gz",
        type=["fastq", "fq", "gz"],
        accept_multiple_files=True,
    )

    st.header("Anchors")
    left_anchor = clean_seq(st.text_input("Left anchor", "CTATAAAAGAGCTCACAACCCCTCA"))
    right_anchor = clean_seq(st.text_input("Right anchor", "GGAGGCCACACCCGCCACTCACCTG"))

    st.header("Parsing options")
    orientation = st.selectbox(
        "Read orientation",
        ["both", "forward_only", "reverse_complement_only"],
        index=0,
    )
    anchor_mismatch = st.number_input("Max mismatch per anchor", min_value=0, max_value=3, value=0, step=1)
    min_len = st.number_input("Minimum insert length", min_value=0, value=0, step=1)
    max_len = st.number_input("Maximum insert length", min_value=1, value=10000, step=1)

    st.header("Optional reference")
    ref_upload = st.file_uploader("Reference CSV", type=["csv"])


ref_path = None
ref_preview = None
id_col = "utr_id"
seq_col = "utr_seq"

if ref_upload is not None:
    ref_preview = pd.read_csv(ref_upload)
    st.subheader("Reference preview")
    st.dataframe(ref_preview.head(), use_container_width=True)
    cols = list(ref_preview.columns)
    col1, col2 = st.columns(2)
    with col1:
        id_col = st.selectbox("Reference ID column", cols, index=0)
    with col2:
        seq_col = st.selectbox("Reference sequence column", cols, index=min(1, len(cols) - 1))

run = st.button("Run analysis", type="primary")

if run:
    if not fastq_files:
        st.error("Upload at least one FASTQ file.")
        st.stop()
    if not left_anchor or not right_anchor:
        st.error("Both left and right anchors are required.")
        st.stop()

    with st.spinner("Running NGS_LibraryQC..."):
        fastq_paths = [save_uploaded_file(f) for f in fastq_files]
        if ref_upload is not None:
            ref_path = save_uploaded_csv(ref_upload)

        insert_df, ref_df, nonref_df, summary = run_counting(
            fastq_paths=fastq_paths,
            left=left_anchor,
            right=right_anchor,
            ref_path=ref_path,
            id_col=id_col,
            seq_col=seq_col,
            orientation=orientation,
            anchor_mismatch=int(anchor_mismatch),
            min_len=int(min_len),
            max_len=int(max_len),
        )

    st.success("Analysis complete")

    st.subheader("Summary")
    summary_df = pd.DataFrame(summary.items(), columns=["metric", "value"])
    st.dataframe(summary_df, use_container_width=True)
    st.download_button(
        "Download summary CSV",
        summary_df.to_csv(index=False).encode(),
        file_name="ngs_libraryqc.summary.csv",
        mime="text/csv",
    )

    tab1, tab2, tab3, tab4 = st.tabs([
        "Extracted inserts",
        "Reference counts",
        "Non-reference",
        "Plots",
    ])

    with tab1:
        st.dataframe(insert_df.head(10000), use_container_width=True)
        st.download_button(
            "Download insert counts",
            insert_df.to_csv(index=False).encode(),
            file_name="ngs_libraryqc.insert_counts.csv",
            mime="text/csv",
        )

    with tab2:
        if ref_df is None:
            st.info("Upload a reference CSV to generate reference-level counts.")
        else:
            st.dataframe(ref_df.sort_values("raw_count", ascending=False), use_container_width=True)
            st.download_button(
                "Download reference counts",
                ref_df.to_csv(index=False).encode(),
                file_name="ngs_libraryqc.reference_counts.csv",
                mime="text/csv",
            )

    with tab3:
        if ref_df is None:
            st.info("Non-reference classification requires a reference CSV.")
        else:
            st.dataframe(nonref_df.head(10000), use_container_width=True)
            st.download_button(
                "Download non-reference sequences",
                nonref_df.to_csv(index=False).encode(),
                file_name="ngs_libraryqc.non_reference.csv",
                mime="text/csv",
            )

    with tab4:
        if ref_df is not None:
            fig, ax = plt.subplots()
            ax.hist(np.log10(ref_df["raw_count"] + 1), bins=50)
            ax.set_xlabel("log10(count + 1)")
            ax.set_ylabel("Number of reference sequences")
            ax.set_title("Reference count distribution")
            st.pyplot(fig)

            ranked = ref_df.sort_values("raw_count", ascending=False).reset_index(drop=True)
            fig2, ax2 = plt.subplots()
            ax2.plot(np.arange(1, len(ranked) + 1), ranked["raw_count"].values)
            ax2.set_xlabel("Rank")
            ax2.set_ylabel("Raw count")
            ax2.set_title("Ranked reference abundance")
            st.pyplot(fig2)
        else:
            fig, ax = plt.subplots()
            ax.hist(np.log10(insert_df["count"] + 1), bins=50)
            ax.set_xlabel("log10(count + 1)")
            ax.set_ylabel("Number of unique inserts")
            ax.set_title("Extracted insert count distribution")
            st.pyplot(fig)
