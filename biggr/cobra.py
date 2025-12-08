import cobra as cobrapy

from biggr import objects
from biggr.models import (
    CompartmentalizedComponent,
    Component,
    UniversalCompartmentalizedComponent,
    UniversalComponent,
)

METABOLITE_ANNOTATION_PRIORITY = ["BiGGr", "BiGG", "CHEBI", "seed.compound"]


def find_metabolite(
    metabolite: cobrapy.Metabolite, cobra_id_namespace="BiGGr", default_compartment=None
):
    model_bigg_id = model.id if (model := metabolite.model) is not None else None

    for ann_type in METABOLITE_ANNOTATION_PRIORITY:
        ann_ids = []
        if ann_type in metabolite.annotation:
            ann_ids.extend(metabolite.annotation[ann_type])
        if ann_type.upper() == cobra_id_namespace.upper():
            ann_ids.append(metabolite.id)
        if not ann_ids:
            continue
        ann_ids = [
            (x if x.upper().startswith(f"{ann_type.upper()}:") else f"{ann_type}:{x}")
            for x in ann_ids
        ]
        m = objects.get_metabolites_by_identifiers(ann_ids, model_bigg_id=model_bigg_id)
        m = [x for x in m.values() if x is not None]
        if not m:
            continue
        result = []

        # Best case is to match a compartmentalized component.
        m_sel = [x for x in m if isinstance(x, CompartmentalizedComponent)]
        for x in m_sel:
            if float(x.component.charge) != float(metabolite.charge):
                continue
            # TODO: Proper formula comparison
            if x.component.formula != metabolite.formula:
                continue
            result.append(x)
        if len(result) == 1:
            return result[0]
        elif len(result) > 1:
            return None

        # Next best case is to match a universal compartmentalized component.
        m_sel = [x for x in m if isinstance(x, UniversalCompartmentalizedComponent)]
        for x in m_sel:
            for cc in x.compartmentalized_components:
                if float(cc.component.charge) != float(metabolite.charge):
                    continue
                # TODO: Proper formula comparison
                if cc.component.formula != metabolite.formula:
                    continue
                result.append(cc)

        if len(result) == 1:
            return result[0]
        elif len(result) > 1:
            return None

        if default_compartment is None:
            return None

        # If a default compartment is specified, we can use a component.
        m_sel = [x for x in m if isinstance(x, Component)]
        for x in m_sel:
            if float(x.charge) != float(metabolite.charge):
                continue
            # TODO: Proper formula comparison
            if x.formula != metabolite.formula:
                continue
            for cc in x.compartmentalized_components:
                if cc.compartment.bigg_id == default_compartment:
                    result.append(cc)

        if len(result) == 1:
            return result[0]
        elif len(result) > 1:
            return None

        # If a default compartment is specified, we can use a universal component.
        m_sel = [x for x in m if isinstance(x, UniversalComponent)]
        for x in m_sel:
            for c in x.components:
                if float(c.charge) != float(metabolite.charge):
                    continue
                # TODO: Proper formula comparison
                if c.formula != metabolite.formula:
                    continue
                for cc in c.compartmentalized_components:
                    if cc.compartment.bigg_id == default_compartment:
                        result.append(cc)

        if len(result) == 1:
            return result[0]
        elif len(result) > 1:
            return None
        return None


def update_metabolite(
    metabolite: cobrapy.Metabolite,
    compartmentalized_component: CompartmentalizedComponent,
):
    bigg_id = compartmentalized_component.bigg_id
    bigg_id = compartmentalized_component.universal_compartmentalized_component.bigg_id
    metabolite.id = bigg_id
    metabolite.name = compartmentalized_component.component.name
    metabolite.annotation.clear()
    metabolite.annotation["BiGGr"] = [compartmentalized_component.bigg_id]

    references = [
        x.reference_compound
        for x in compartmentalized_component.component.reference_mappings
    ]
    for reference in references:
        # if reference.bigg_id.startswith("CHEBI:"):
        #     if "CHEBI" in metabolite.annotation:
        #         metabolite.annotation["CHEBI"].append(reference.bigg_id)
        #     else:
        #         metabolite.annotation["CHEBI"] = [reference.bigg_id]
        if reference.inchi_id is not None:
            # if "InChI" in metabolite.annotation:
            #     metabolite.annotation["InChI"].append(reference.inchi.to_string())
            # else:
            #     metabolite.annotation["InChI"] = reference.inchi.to_string()
            if "InChIKey" in metabolite.annotation:
                metabolite.annotation["InChIKey"].append(reference.inchi.key())
            else:
                metabolite.annotation["InChIKey"] = [reference.inchi.key()]
        for ann_map in reference.annotation_mappings:
            annotation = ann_map.annotation
            for link in annotation.links:
                ds = link.data_source
                namespace = ds.bigg_id
                identifier = link.identifier
                if namespace in metabolite.annotation:
                    if identifier not in metabolite.annotation[namespace]:
                        metabolite.annotation[namespace].append(identifier)
                else:
                    metabolite.annotation[namespace] = [identifier]
    for ann_map in compartmentalized_component.component.annotation_mappings:
        annotation = ann_map.annotation
        for link in annotation.links:
            ds = link.data_source
            namespace = ds.bigg_id
            identifier = link.identifier
            if namespace in metabolite.annotation:
                if identifier not in metabolite.annotation[namespace]:
                    metabolite.annotation[namespace].append(identifier)
            else:
                metabolite.annotation[namespace] = [identifier]

    return metabolite
