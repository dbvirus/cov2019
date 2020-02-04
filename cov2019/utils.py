"""
tool box of useful functions
"""
import logging
import os

from dbvirus_searcher.searcher import Searcher


def srr_to_url(srr, halt_on_failure=False):
    """
    Given an SRA SRR number, returns its URL.

    In case of failure, it returns None.

    Please provide an e-mail address in the NCBI_EMAIL environment variable. It
    is not necessary, but it is a good practice.

    If `halt_on_failure` is set to True, an Exception is thrown if
    any errors are found while acquiring data.

    Example of usage:
    >>> srr_to_url("SRR042301")
    'https://sra-download.st-va.ncbi.nlm.nih.gov/sos2/sra-pub-run-6/SRR042301/SRR042301.3'

    >>> srr_to_url("invalid")
    
    >>> srr_to_url("042301")
    'https://sra-download.st-va.ncbi.nlm.nih.gov/sos2/sra-pub-run-6/SRR042301/SRR042301.3'

    >>> srr_to_url("SRR327703")
    'https://sra-download.st-va.ncbi.nlm.nih.gov/sos2/sra-pub-run-3/SRR327703/SRR327703.3'

    >>> srr_to_url("SRR10059476")
    'https://sra-download.ncbi.nlm.nih.gov/traces/sra66/SRR/009823/SRR10059476'

    >>> srr_to_url("SRR10059477")
    'https://sra-download.ncbi.nlm.nih.gov/traces/sra35/SRR/009823/SRR10059477'
    """

    srr = srr.upper()

    email = os.environ.get("NCBI_EMAIL", "default@ncbi.com")

    if email == "default@ncbi.com":
        logging.warn(
            "Please provide a custom email address by setting the NCBI_EMAIL env"
        )

    # If the SRR string does not start with SRR, we add it
    if not srr.startswith("SRR"):
        srr = f"SRR{srr}"

    # Here we perform the basic SRA search
    searcher = Searcher(email=email)
    results = searcher.search(query=srr).get("esearchresult")
    id_list = results.get("idlist")

    if not id_list:
        return None

    id_ = id_list[0]
    data = searcher.fetch(id_)

    # We should not have more than 1 UID. Ever. Just in case we do, log it!
    if len(id_list) > 1:
        logging.warn(f"More than 1 UID found for {srr}")
        logging.warn(f"Processing only the first UID, {id_}")
    
    # The result has a very complicated nested structure.
    # Here we try to access the target data.
    try:
        runs = data["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]["RUN_SET"]["RUN"]

        # It is possible that we are dealing with a list of runs
        if isinstance(runs, list):
            # If so, there is only one accession # that matches the target SRR
            accessions = [r["@accession"] for r in runs]
            run_id = accessions.index(srr)
            run = runs[run_id]
        else:
            # If not, there's a single run
            run = runs
        
        sra_file = run["SRAFiles"]["SRAFile"]

        # A single SRA File may have multiple REAL files
        # See SRR 10059477 for a current example
        if isinstance(sra_file, list):
            # Currently, we care about only for those with semantic_name = run
            semantic_names = [f["@semantic_name"] for f in sra_file]
            file_id = semantic_names.index("run")

            if not file_id:
                logging.warn(f"SRR {srr} has multiple files but no one is valid!")
                logging.warn("No file has 'semantic_name = run' ")
                return None

            sra_file = sra_file[file_id]

        url = sra_file["@url"]
    except Exception as e:
        logging.error(f"Failed to acquire link for SRR {srr}")
        logging.error(e)

        if halt_on_failure:
            raise Exception(f"Failed to acquire link for SRR {srr}")

        return None
    else:
        return url

