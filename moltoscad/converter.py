"""MOL to SCAD conversion utilities."""

import os
import subprocess
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem


__all__ = [
    "add_hydrogens",
    "mol_to_scad",
    "get_atom_color",
    "get_atom_radius",
    "open_scad_file",
    "open_freecad_with_scad",
]


def add_hydrogens(mol: Chem.Mol) -> Chem.Mol:
    """Add explicit hydrogens and generate 3D coordinates."""
    mol = Chem.AddHs(mol)

    methods = [
        (AllChem.EmbedMolecule, {"randomSeed": 42, "useRandomCoords": True}),
        (AllChem.EmbedMolecule, {"randomSeed": 42, "useRandomCoords": True, "useBasicKnowledge": False}),
        (AllChem.EmbedMolecule, {"randomSeed": 42, "useRandomCoords": True, "boxSizeMult": 2.0}),
        (AllChem.EmbedMolecule, {"randomSeed": 42, "useExpTorsionAnglePrefs": True, "useBasicKnowledge": True}),
    ]

    for method, params in methods:
        result = method(mol, **params)
        if result == 0:
            AllChem.MMFFOptimizeMolecule(mol)
            return mol

    raise ValueError("Could not generate 3D coordinates")


def mol_to_scad(mol_file: str, scad_file: str, max_atoms_per_file: int = 1000) -> bool:
    """Convert a MOL file to a SCAD file."""
    with open(mol_file, "r") as f:
        mol_block = f.read()
    mol = Chem.MolFromMolBlock(mol_block, removeHs=False)

    if mol is None:
        raise ValueError(f"Could not read MOL file: {mol_file}")

    mol = Chem.RemoveHs(mol)
    mol = Chem.AddHs(mol)
    mol = add_hydrogens(mol)

    conf = mol.GetConformer()
    if not conf.Is3D():
        raise ValueError("3D coordinates were not generated")

    atom_count = mol.GetNumAtoms()
    file_count = (atom_count - 1) // max_atoms_per_file + 1

    for file_index in range(file_count):
        start_atom = file_index * max_atoms_per_file
        end_atom = min((file_index + 1) * max_atoms_per_file, atom_count)

        current_scad_file = (
            f"{os.path.splitext(scad_file)[0]}_{file_index + 1}.scad"
            if file_count > 1
            else scad_file
        )

        with open(current_scad_file, "w") as f:
            f.write("// Generated SCAD file from MOL\n\n")

            f.write("// Adjustable Parameters\n")
            f.write("atom_scale = 1;  // Adjusting the size of atoms\n")
            f.write("bond_radius = 0.2;  // Coupling Radius\n\n")

            f.write("// Atomic color\n")
            for element, color in get_atom_color.colors.items():
                f.write(f"{element}_color = \"{color}\";\n")
            f.write("\n")

            f.write("// Atomic radius\n")
            for element, radius in get_atom_radius.radii.items():
                f.write(f"{element}_radius = {radius};\n")
            f.write("\n")

            for atom_idx in range(start_atom, end_atom):
                atom = mol.GetAtomWithIdx(atom_idx)
                pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
                element = atom.GetSymbol()

                f.write(
                    f"translate([{pos.x:.2f}, {pos.y:.2f}, {pos.z:.2f}]) "
                    f"color({element}_color) sphere(r = {element}_radius * atom_scale);\n"
                )

            for bond in mol.GetBonds():
                atom1_idx = bond.GetBeginAtomIdx()
                atom2_idx = bond.GetEndAtomIdx()
                if start_atom <= atom1_idx < end_atom or start_atom <= atom2_idx < end_atom:
                    pos1 = mol.GetConformer().GetAtomPosition(atom1_idx)
                    pos2 = mol.GetConformer().GetAtomPosition(atom2_idx)
                    f.write(
                        "color(\"grey\") "
                        f"translate([{pos1.x:.2f}, {pos1.y:.2f}, {pos1.z:.2f}]) "
                        f"rotate([0,acos({pos2.z-pos1.z:.2f}/norm([{pos2.x-pos1.x:.2f},{pos2.y-pos1.y:.2f},{pos2.z-pos1.z:.2f}])),atan2({pos2.y-pos1.y:.2f},{pos2.x-pos1.x:.2f})]) "
                        f"cylinder(h = norm([{pos2.x-pos1.x:.2f},{pos2.y-pos1.y:.2f},{pos2.z-pos1.z:.2f}]), r = bond_radius);\n"
                    )

    return True


def get_atom_color(element: str) -> str:
    return get_atom_color.colors.get(element, "pink")


get_atom_color.colors = {
    "H": "white",
    "C": "black",
    "N": "navy",
    "O": "red",
    "F": "green",
    "Cl": "green",
    "Br": "brown",
    "I": "darkviolet",
    "He": "cyan",
    "Ne": "cyan",
    "Ar": "cyan",
    "Xe": "cyan",
    "Kr": "cyan",
    "P": "orange",
    "S": "yellow",
    "B": "pink",
    "Li": "purple",
    "Na": "purple",
    "K": "purple",
    "Rb": "purple",
    "Cs": "purple",
    "Be": "darkgreen",
    "Mg": "darkgreen",
    "Ca": "darkgreen",
    "Sr": "darkgreen",
    "Ba": "darkgreen",
    "Ra": "darkgreen",
    "Ti": "gray",
    "Fe": "orange",
}


def get_atom_radius(element: str) -> float:
    original = get_atom_radius.original_radii.get(element, 0.77)
    return original * get_atom_radius.normalization_factor


get_atom_radius.original_radii = {
    "H": 0.37,
    "He": 0.32,
    "Li": 1.34,
    "Be": 0.9,
    "B": 0.82,
    "C": 0.77,
    "N": 0.75,
    "O": 0.73,
    "F": 0.71,
    "Ne": 0.69,
    "Na": 1.54,
    "Mg": 1.3,
    "Al": 1.18,
    "Si": 1.11,
    "P": 1.06,
    "S": 1.02,
    "Cl": 0.99,
    "Ar": 0.97,
    "K": 1.96,
    "Ca": 1.74,
    "Sc": 1.44,
    "Ti": 1.36,
    "V": 1.25,
    "Cr": 1.27,
    "Mn": 1.39,
    "Fe": 1.25,
    "Co": 1.26,
    "Ni": 1.21,
    "Cu": 1.38,
    "Zn": 1.31,
    "Ga": 1.26,
    "Ge": 1.22,
    "As": 1.19,
    "Se": 1.16,
    "Br": 1.14,
    "Kr": 1.1,
    "Rb": 2.11,
    "Sr": 1.92,
    "Y": 1.62,
    "Zr": 1.48,
    "Nb": 1.37,
    "Mo": 1.45,
    "Tc": 1.56,
    "Ru": 1.26,
    "Rh": 1.35,
    "Pd": 1.31,
    "Ag": 1.53,
    "Cd": 1.48,
    "In": 1.44,
    "Sn": 1.41,
    "Sb": 1.38,
    "Te": 1.35,
    "I": 1.33,
    "Xe": 1.3,
    "Cs": 2.11,
    "Ba": 1.92,
    "Ce": 1.62,
    "Pr": 1.48,
    "Nd": 1.37,
    "Pm": 1.45,
    "Sm": 1.56,
    "Eu": 1.26,
    "Gd": 1.35,
    "Tb": 1.31,
    "Dy": 1.53,
    "Ho": 1.48,
    "Er": 1.44,
    "Tm": 1.41,
    "Yb": 1.38,
    "Lu": 1.35,
    "Hf": 1.33,
    "Ta": 1.3,
    "W": 1.3,
    "Re": 1.3,
    "Os": 1.3,
    "Ir": 1.3,
    "Pt": 1.3,
    "Au": 1.3,
    "Hg": 1.3,
    "Tl": 1.3,
    "Pb": 1.3,
    "Bi": 1.3,
    "Po": 1.3,
    "At": 1.3,
    "Rn": 1.3,
    "Fr": 2.11,
    "Ra": 1.92,
    "Ac": 1.62,
    "Pa": 1.48,
    "U": 1.37,
    "Np": 1.45,
    "Pu": 1.56,
    "Am": 1.26,
    "Cm": 1.35,
    "Bk": 1.31,
    "Cf": 1.53,
    "Es": 1.48,
    "Fm": 1.44,
    "Md": 1.41,
    "No": 1.38,
    "Lr": 1.35,
}


carbon_radius = get_atom_radius.original_radii["C"]
get_atom_radius.normalization_factor = 0.5 / carbon_radius
get_atom_radius.radii = {
    element: get_atom_radius(element)
    for element in get_atom_radius.original_radii
}


def open_scad_file(scad_file: str, openscad_path: Optional[str] = None) -> None:
    if openscad_path is None:
        openscad_path = r"C:\\Program Files (x86)\\OpenSCAD\\openscad.exe"
    try:
        subprocess.Popen([openscad_path, scad_file])
    except FileNotFoundError:
        raise FileNotFoundError("OpenSCAD not found. Check the path setting.")


def open_freecad_with_scad(scad_file: str, freecad_path: Optional[str] = None) -> None:
    if freecad_path is None:
        freecad_path = r"C:\\Program Files\\FreeCAD 0.21\\bin\\FreeCAD.exe"
    try:
        subprocess.Popen([freecad_path, scad_file])
    except FileNotFoundError:
        raise FileNotFoundError("FreeCAD not found. Check the path setting.")
