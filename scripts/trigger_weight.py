from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import uproot
from correctionlib import schemav2

SCRIPT_PATH = Path(__file__).parent
DATA_PATH = SCRIPT_PATH / "data"
DATA_PATH.mkdir(parents=True, exist_ok=True)


def get_corr_pt_msd(eff, eff_unc_up, eff_unc_dn, year, label, trigger_label, edges):
    def multibinning(eff):
        return schemav2.MultiBinning(
            nodetype="multibinning",
            inputs=["pt", "msd"],
            edges=edges,
            content=list(eff.flatten()),
            flow=1.0,  # SET FLOW TO 1.0
        )

    corr = schemav2.Correction(
        name=f"fatjet_triggereff{label}_{year}_{trigger_label}",
        description=f"{label} Trigger efficiency for trigger soup: {trigger_label}",
        version=1,
        inputs=[
            schemav2.Variable(
                name="pt",
                type="real",
                description="Jet transverse momentum",
            ),
            schemav2.Variable(
                name="msd",
                type="real",
                description="Jet softdrop mass",
            ),
            schemav2.Variable(
                name="systematic",
                type="string",
                description="Systematic variation",
            ),
        ],
        output=schemav2.Variable(
            name="weight", type="real", description=f"Jet {label} trigger efficiency"
        ),
        data=schemav2.Category(
            nodetype="category",
            input="systematic",
            content=[
                {"key": "nominal", "value": multibinning(eff)},
                {"key": "stat_up", "value": multibinning(eff_unc_up)},
                {"key": "stat_dn", "value": multibinning(eff_unc_dn)},
            ],
        ),
    )
    return corr


def get_corr_txbb(eff, eff_unc_up, eff_unc_dn, year, label, trigger_label, edges):
    def singlebinning(eff):
        return schemav2.Binning(
            nodetype="binning",
            input="txbb",
            edges=edges,
            content=list(eff.flatten()),
            flow=1.0,  # SET FLOW TO 1.0
        )

    corr = schemav2.Correction(
        name=f"fatjet_triggereff{label}_{year}_{trigger_label}",
        description=f"{label} Trigger efficiency for trigger soup: {trigger_label}",
        version=1,
        inputs=[
            schemav2.Variable(
                name="txbb",
                type="real",
                description="Jet TXbb",
            ),
            schemav2.Variable(
                name="systematic",
                type="string",
                description="Systematic variation",
            ),
        ],
        output=schemav2.Variable(
            name="weight", type="real", description=f"Jet {label} trigger efficiency"
        ),
        data=schemav2.Category(
            nodetype="category",
            input="systematic",
            content=[
                {"key": "nominal", "value": singlebinning(eff)},
                {"key": "stat_up", "value": singlebinning(eff_unc_up)},
                {"key": "stat_dn", "value": singlebinning(eff_unc_dn)},
            ],
        ),
    )
    return corr


def main(args):
    region = args.region

    corr_pt_msd = {
        "2022": {
            "fname": f"sf_mass_pt/output/efficiency_mass_pt_2022_{region}.root",
            "data": "eff_data",
            "mc": "eff_MC",
        },
        "2022EE": {
            "fname": f"sf_mass_pt/output/efficiency_mass_pt_2022EE_{region}.root",
            "data": "eff_data",
            "mc": "eff_MC",
        },
        "2023": {
            "fname": f"sf_mass_pt/output/efficiency_mass_pt_2023_{region}.root",
            "data": "eff_data",
            "mc": "eff_MC",
        },
        "2023BPix": {
            "fname": f"sf_mass_pt/output/efficiency_mass_pt_2023BPix_{region}.root",
            "data": "eff_data",
            "mc": "eff_MC",
        },
    }

    corr_txbb_glopart = {
        "2022": {
            "fname": f"sf_TXbb/output/efficiency_GloParT_2022_{region}.root",
            "data": "eff_data",
            "mc": "eff_MC",
        },
        "2022EE": {
            "fname": f"sf_TXbb/output/efficiency_GloParT_2022EE_{region}.root",
            "data": "eff_data",
            "mc": "eff_MC",
        },
        "2023": {
            "fname": f"sf_TXbb/output/efficiency_GloParT_2023_{region}.root",
            "data": "eff_data",
            "mc": "eff_MC",
        },
        "2023BPix": {
            "fname": f"sf_TXbb/output/efficiency_GloParT_2023BPix_{region}.root",
            "data": "eff_data",
            "mc": "eff_MC",
        },
    }

    corr_txbb_pnet = {
        "2022": {
            "fname": f"sf_TXbb/output/efficiency_PNetLegacy_2022_{region}.root",
            "data": "eff_data",
            "mc": "eff_MC",
        },
        "2022EE": {
            "fname": f"sf_TXbb/output/efficiency_PNetLegacy_2022EE_{region}.root",
            "data": "eff_data",
            "mc": "eff_MC",
        },
        "2023": {
            "fname": f"sf_TXbb/output/efficiency_PNetLegacy_2023_{region}.root",
            "data": "eff_data",
            "mc": "eff_MC",
        },
        "2023BPix": {
            "fname": f"sf_TXbb/output/efficiency_PNetLegacy_2023BPix_{region}.root",
            "data": "eff_data",
            "mc": "eff_MC",
        },
    }

    for year, corr_items in corr_pt_msd.items():
        f = uproot.open(corr_items["fname"])
        corrs = {}
        for key in ["data", "mc"]:
            h = f[corr_items[key]].to_hist()
            # Since input is (mass, pt), we need to transpose and swap edges
            edges = [
                h.axes[1].edges,  # pt edges
                h.axes[0].edges,  # mass edges
            ]

            # Transpose the values to match the desired (pt, mass) order
            values = h.values().T
            variances = h.variances().T

            corr = get_corr_pt_msd(
                values,
                np.sqrt(variances),
                np.sqrt(variances),
                year,
                key,
                "ptmsd",
                edges,
            )
            corrs[key] = corr

        cset = schemav2.CorrectionSet(schema_version=2, corrections=[corrs["mc"], corrs["data"]])
        path = Path(DATA_PATH / f"fatjet_triggereff_{year}_ptmsd_{region}.json")
        with path.open("w") as fout:
            fout.write(cset.json(exclude_unset=True))

    for year, corr_items in corr_txbb_glopart.items():
        f = uproot.open(corr_items["fname"])
        corrs = {}
        for key in ["data", "mc"]:
            h = f[corr_items[key]].to_hist()
            edges = h.axes.edges[0]
            corr = get_corr_txbb(
                h.values(),
                np.sqrt(h.variances()),
                np.sqrt(h.variances()),
                year,
                key,
                "txbbGloParT",
                edges,
            )
            corrs[key] = corr

        cset = schemav2.CorrectionSet(schema_version=2, corrections=[corrs["mc"], corrs["data"]])
        path = Path(DATA_PATH / f"fatjet_triggereff_{year}_txbbGloParT_{region}.json")
        with path.open("w") as fout:
            fout.write(cset.json(exclude_unset=True))

    for year, corr_items in corr_txbb_pnet.items():
        f = uproot.open(corr_items["fname"])
        corrs = {}
        for key in ["data", "mc"]:
            h = f[corr_items[key]].to_hist()
            edges = h.axes.edges[0]
            corr = get_corr_txbb(
                h.values(),
                np.sqrt(h.variances()),
                np.sqrt(h.variances()),
                year,
                key,
                "txbbPNetLegacy",
                edges,
            )
            corrs[key] = corr

        cset = schemav2.CorrectionSet(schema_version=2, corrections=[corrs["mc"], corrs["data"]])
        path = Path(DATA_PATH / f"fatjet_triggereff_{year}_txbbPNet_{region}.json")
        with path.open("w") as fout:
            fout.write(cset.json(exclude_unset=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--region", required=True, choices=["QCD", "TTbar"], type=str)
    args = parser.parse_args()

    main(args)