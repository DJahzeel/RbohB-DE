
from Bio import Entrez

def explore_databases(project_uid, database):
    handle = Entrez.elink(dbfrom="bioproject", db=database, id=project_uid)
    linksBioProj = Entrez.read(handle)
    handle.close()
    return linksBioProj

def geo(linksBioProj,project_accession):
    if linksBioProj[0]["LinkSetDb"]:
        geo_ids = [link["Id"] for link in linksBioProj[0]["LinkSetDb"][0]["Link"]]
        print("\nEl proyecto", project_accession, " está asociado a la base de datos de GEO")
        print("UIDs GEO asociados:", geo_ids[:5], "... total:", len(geo_ids))
        for geo_id in geo_ids:
            handle = Entrez.esummary(db="gds", id=geo_id)
            summary = Entrez.read(handle)
            handle.close()
            for key, value in summary[0].items():
                print(f"{key}: {value}")
    else:
        return None
    return geo_ids

def sra(linksBioProj,project_accession):
    if linksBioProj and linksBioProj [0].get("LinkSetDb") and linksBioProj[0]["LinkSetDb"]:
        sra_ids = [link["Id"] for link in linksBioProj[0]["LinkSetDb"][0]["Link"]]
        print("\nEl proyecto", project_accession, " está asociado a la base de datos de SRA")
        print("UIDs SRA asociados:", sra_ids, "... total:", len(sra_ids))
    else:
        return None
    return sra_ids

