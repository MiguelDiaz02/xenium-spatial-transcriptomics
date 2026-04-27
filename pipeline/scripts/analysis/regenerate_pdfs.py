#!/usr/bin/env python3
"""
Regenerate all Phase 3 PDFs with updated figures
"""

import time
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from PIL import Image
import logging

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')
log = logging.getLogger(__name__)

FIG_DIR_ENRICH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_neighborhood_enrichment"
FIG_DIR_REGION_REF = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_tissue_region_refinement"


def regenerate_neighborhood_enrichment_pdf():
    """Regenerate Phase3_Task3.2_Neighborhood_Enrichment.pdf"""
    log.info("Regenerating Phase 3.2 PDF (Neighborhood Enrichment)...")

    pdf_path = FIG_DIR_ENRICH / "Phase3_Task3.2_Neighborhood_Enrichment.pdf"

    with PdfPages(pdf_path) as pdf:
        # Title page
        fig = plt.figure(figsize=(11, 8.5))
        ax = fig.add_subplot(111)
        ax.axis('off')
        ax.text(0.5, 0.7, 'Phase 3.2: Neighborhood Enrichment', fontsize=24, fontweight='bold',
               ha='center', transform=ax.transAxes)
        ax.text(0.5, 0.6, 'Cell Type Co-localization Patterns', fontsize=16,
               ha='center', transform=ax.transAxes)
        ax.text(0.5, 0.5, '268,034 cells × 7 cell types', fontsize=12,
               ha='center', transform=ax.transAxes, style='italic')
        ax.text(0.5, 0.3, 'Updated with corrected figures', fontsize=10,
               ha='center', transform=ax.transAxes, color='green', fontweight='bold')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Figure 2: Enrichment Patterns (UPDATED - taller)
        img_path = FIG_DIR_ENRICH / "Fig2_Enrichment_Patterns.png"
        if img_path.exists():
            img = Image.open(img_path)
            fig = plt.figure(figsize=(11, 14))
            ax = fig.add_subplot(111)
            ax.imshow(img)
            ax.axis('off')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Figure 3: Spatial Localization (UPDATED - populated)
        img_path = FIG_DIR_ENRICH / "Fig3_Spatial_Localization.png"
        if img_path.exists():
            img = Image.open(img_path)
            fig = plt.figure(figsize=(11, 9))
            ax = fig.add_subplot(111)
            ax.imshow(img)
            ax.axis('off')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Figure 4: Spatial Niches (UPDATED - tiny points)
        img_path = FIG_DIR_ENRICH / "Fig4_Spatial_Niches.png"
        if img_path.exists():
            img = Image.open(img_path)
            fig = plt.figure(figsize=(11, 8))
            ax = fig.add_subplot(111)
            ax.imshow(img)
            ax.axis('off')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Summary page
        fig = plt.figure(figsize=(11, 8.5))
        ax = fig.add_subplot(111)
        ax.axis('off')
        summary_text = """
        ANALYSIS SUMMARY

        ✓ Enrichment Patterns: Complete 7×7 cell type matrix
          - Shows Log Odds Ratio of spatial co-localization
          - Identifies cell type clustering patterns

        ✓ Spatial Localization: Pairwise spatial patterns
          - T_cell vs Macrophage, T_cell vs Tumor, etc.
          - Validates co-localization enrichment

        ✓ Spatial Niches: Functional compartments
          - Microscale organization of immune microenvironment
          - Links enrichment to spatial organization

        All figures updated: 2026-04-27
        """
        ax.text(0.1, 0.9, summary_text, fontsize=10, verticalalignment='top',
               transform=ax.transAxes, family='monospace',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    log.info(f"  ✓ Saved: {pdf_path}")


def regenerate_boundary_trajectories_pdf():
    """Regenerate Phase3_Task3.4_Tissue_Region_Refinement.pdf"""
    log.info("Regenerating Phase 3.4 PDF (Boundary Validation & Trajectories)...")

    pdf_path = FIG_DIR_REGION_REF / "Phase3_Task3.4_Tissue_Region_Refinement.pdf"

    with PdfPages(pdf_path) as pdf:
        # Title page
        fig = plt.figure(figsize=(11, 8.5))
        ax = fig.add_subplot(111)
        ax.axis('off')
        ax.text(0.5, 0.7, 'Phase 3.4: Tissue Region Refinement', fontsize=24, fontweight='bold',
               ha='center', transform=ax.transAxes)
        ax.text(0.5, 0.6, 'Advanced Analysis: Boundaries & Transitions', fontsize=16,
               ha='center', transform=ax.transAxes)
        ax.text(0.5, 0.5, '3 tissue regions × 6 advanced analyses', fontsize=12,
               ha='center', transform=ax.transAxes, style='italic')
        ax.text(0.5, 0.3, 'Updated with corrected figures', fontsize=10,
               ha='center', transform=ax.transAxes, color='green', fontweight='bold')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Figure 4: Boundary Validation (UPDATED - legend outside)
        img_path = FIG_DIR_REGION_REF / "Fig4_Boundary_Validation.png"
        if img_path.exists():
            img = Image.open(img_path)
            fig = plt.figure(figsize=(11, 5))
            ax = fig.add_subplot(111)
            ax.imshow(img)
            ax.axis('off')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Figure 5: Trajectories (UPDATED - spacing + legend)
        img_path = FIG_DIR_REGION_REF / "Fig5_Spatial_Trajectories.png"
        if img_path.exists():
            img = Image.open(img_path)
            fig = plt.figure(figsize=(11, 5))
            ax = fig.add_subplot(111)
            ax.imshow(img)
            ax.axis('off')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Summary page
        fig = plt.figure(figsize=(11, 8.5))
        ax = fig.add_subplot(111)
        ax.axis('off')
        summary_text = """
        REGION REFINEMENT SUMMARY

        ✓ Boundary Validation (3 regions):
          - Internal cohesion: High within-region consistency
          - Boundary clarity: Sharp region definitions
          - Spatial dispersion: Characteristic region sizes

        ✓ Spatial Trajectories (3 pairs):
          - Expression correlation between adjacent regions
          - Transition sharpness: Smooth vs abrupt boundaries
          - Links region organization to molecular patterns

        ✓ Additional Analyses:
          - Pathway enrichment (10 pathways × 3 regions)
          - Immune signature scoring (7 signatures)
          - DE ranking (all 289 genes per region)

        All figures updated: 2026-04-27
        """
        ax.text(0.1, 0.9, summary_text, fontsize=10, verticalalignment='top',
               transform=ax.transAxes, family='monospace',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    log.info(f"  ✓ Saved: {pdf_path}")


def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("REGENERATING ALL PHASE 3 PDFs WITH UPDATED FIGURES")
    log.info("=" * 80)

    try:
        regenerate_neighborhood_enrichment_pdf()
        regenerate_boundary_trajectories_pdf()

        elapsed = time.time() - start_time
        log.info("\n" + "=" * 80)
        log.info(f"✓ ALL PDFs REGENERATED in {elapsed:.1f} seconds")
        log.info("=" * 80)

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
