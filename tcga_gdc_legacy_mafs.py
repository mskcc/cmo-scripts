#!/usr/bin/env python

import os, requests, json, argparse

def main():
    # The GDC Legacy Archive endpoints we'll need here
    proj_ep = 'https://gdc-api.nci.nih.gov/legacy/projects/'
    file_ep = 'https://gdc-api.nci.nih.gov/legacy/files/'

    parser = argparse.ArgumentParser( description='Download open access TCGA MAFs from GDC Legacy Archive' )
    parser.add_argument('-d', '--download', action='store_true', default=False, help='Specify this if you want to download the files too' )
    args = vars( parser.parse_args() )

    # Fetch all projects IDs, and limit them to anything starting with "TCGA-"
    response = requests.post( proj_ep, data={ 'size':'1000', 'fields':'project_id' })
    proj_ids = [ d['project_id'] for d in response.json()['data']['hits']]
    proj_ids = [ d for d in proj_ids if d.startswith('TCGA-')]

    # For each TCGA project, construct and POST a request that finds all its MAFs
    for proj_id in proj_ids:
        payload = {
            'filters':{
                'op':'and',
                'content':[
                    {
                        'op':'=',
                        'content':{
                            'field':'cases.project.project_id',
                            'value':proj_id
                        }
                    },
                    {
                        'op':'=',
                        'content':{
                            'field':'files.data_format',
                            'value':'MAF'
                        }
                    },
                    {
                        'op':'=',
                        'content':{
                            'field':'access',
                            'value':'open'
                        }
                    }
                ]
            },
            'format':'json',
            'fields':'file_id,file_name,file_size,md5sum',
            'size':'1000'
        }

        response = requests.post( file_ep, json=payload )
        for maf in response.json()['data']['hits']:
            # Print information on the MAFs, and if explicitly requested, download them too
            print( '\t'.join([ proj_id, maf['file_id'], maf['file_name'] ]))
            if args['download']:
               download_file( maf['file_id'], '/'.join([ 'mafs', proj_id, maf['file_name'] ]))

        exit( 0 )

# Given a GDC UUID and a local file path, this function downloads and saves a file
def download_file( file_id, filepath ):
    # Start a stream to the GDC Legacy Archive data endpoint
    data_ep = 'https://gdc-api.nci.nih.gov/legacy/data/'
    response = requests.get( data_ep + file_id, stream=True )

    # Make the folder containing the local file if necessary
    directory = os.path.dirname( filepath )
    if not os.path.exists( directory ):
        os.makedirs( directory )
    # Read and write the file in chunks, which works well for ginormous files
    with open( filepath, 'wb' ) as fh:
        for chunk in response.iter_content( chunk_size=1024 ):
            if chunk: # filter out keep-alive new chunks
                fh.write( chunk )
    response.close()
    return filepath

if __name__ == "__main__":
    main()
