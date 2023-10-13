#!/usr/bin/python

from Bio import Entrez
import xml.etree.ElementTree as ET

Entrez.email = "your.email@mail.com"
Entrez.api_key="entrezkey, looks like hash"

def get_ids(term, db):
    ids = []
    handle = Entrez.esearch(db=db, term=term)
    record = Entrez.read(handle)
    ids.append(record["IdList"])
    return ids

def get_assembly_summary(id):
    handle = Entrez.esummary(db="assembly",id=id,report="full")
    record = Entrez.read(handle)
    #Convert raw output to json
    return(record)

def get_nuccore_summary(id):
    handle = Entrez.efetch(db="nuccore",id=id,retmode="xml")
    record = Entrez.read(handle)
    #Convert raw output to json
    return(record)


def get_SRA_summary(id):
    handle = Entrez.esummary(db="sra",id=id,report="full")
    record = Entrez.read(handle)
    #Convert raw output to json
    return(record)

def get_biosample_summary(id):
    handle = Entrez.esearch(db="biosample", term=id)
    record = Entrez.read(handle)
    
    handle = Entrez.esummary(db="biosample",id=record["IdList"][0], retmode="xml", rettype="full")
    record = Entrez.read(handle)


    bio = record["DocumentSummarySet"]["DocumentSummary"][0]["SampleData"]
    doc = ET.fromstring(bio)
    geoLoc=""
    colDate=""
    for item in doc.findall('.//Attributes//Attribute'):
        for attribute in item.attrib:
            if "harmonized_name" in item.attrib:
                if item.attrib["harmonized_name"]=="geo_loc_name":
                    geoLoc=item.text
                if item.attrib["harmonized_name"]=="collection_date":
                    colDate=item.text
    return [geoLoc, colDate]

#Test
validIDs=[]
validIDs=["SRR9067224"] #this the list of IDs to procecess
startDB="sra" #assembly, sra, nuccore. The database in which to search the above validIDs
output_data="assembly" # or assembly. Sometimes, assembly is not linked to SRA, so in case assembly is not found, nuccore will be output if found

counter=0
with open("~/error_file.tsv","w") as errorsFile:
    with open("~/results.tsv","w") as outputFile:
        for searchterm in validIDs:
            if counter % 10==0:
                print(counter/len(validIDs))
            counter+=1
            try:
                #searchterm="SRR11100134"
                for id in get_ids(searchterm,startDB): #assembly, sra, nuccore, the id is "uid" which can be searched against any database
                    if searchterm=="Run":
                        continue
                    if startDB=="sra":
                        assemblyMetadata=get_SRA_summary(id) #JSON Formatted
                        doc=ET.fromstring("<!DOCTYPE html><html>"+str(assemblyMetadata[0]["ExpXml"])+"</html>")#
                        if output_data=="biosample":
                            output_data=get_biosample_summary( doc.findall('Sample')[0].attrib["acc"] )
                        elif output_data=="assembly": 
                            try: 
                                output_data=get_assembly_summary( doc.findall('Sample')[0].attrib["acc"] )
                            except:
                                output_data=[f["GBSeq_accession-version"] for f in get_nuccore_summary(id)] #there often will be multple IDs

                        #GBSeq_accession-version
                    elif startDB=="assembly":
                        assemblyMetadata=get_assembly_summary(id) #JSON Formatted
                        output_data=get_biosample_summary(assemblyMetadata["DocumentSummarySet"]["DocumentSummary"][0]["BioSampleAccn"])
                    elif startDB=="nuccore":
                        assemblyMetadata=get_nuccore_summary(id) #JSON Formatted
                        output_data=["Unknown","Unknown","Unknown"]
                        for value in assemblyMetadata[0]["GBSeq_feature-table"][0]["GBFeature_quals"]:
                            if value["GBQualifier_name"]=="country":
                                output_data[0]=value["GBQualifier_value"]
                            if value["GBQualifier_name"]=="collection_date":
                                output_data[1]=value["GBQualifier_value"]
                            if value["GBQualifier_name"]=="organism":
                                output_data[2]=value["GBQualifier_value"]
                        print(output_data)
                    outputFile.write('\t'.join( [str(f) for f in [searchterm]+output_data])+"\n")
            except:
                errorsFile.write(searchterm+"\n")

            