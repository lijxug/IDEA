#!/usr/bin/env python3

import sys
import requests
import json
import time
import urllib
import click


class EnsemblRestClient(object):
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urllib.parse.urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        
        try:
	        request = requests.get(self.server + endpoint, headers=hdrs)
	        # if not request.ok:
	        # 	request.raise_for_status()
	        # 	# sys.exit()

	        decoded = request.json()
	        if decoded:
	            data = decoded
	        self.req_count += 1

        except urllib.error.HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))
           
        return data

    def get_variants(self, species, symbol):
        genes = self.perform_rest_action(
            '/xrefs/symbol/{0}/{1}'.format(species, symbol), 
            params={'object_type': 'gene'}
        )
        if genes:
            stable_id = genes[0]['id']
            variants = self.perform_rest_action(
                '/overlap/id/{0}'.format(stable_id),
                params={'feature': 'variation'}
            )
            return variants
        return None

    def get_sequence(self, seq_id, seq_type):
    	sequence = self.perform_rest_action(
    		'/sequence/id/{0}'.format(seq_id),
            params={'type': seq_type}
    	)
    	return sequence

    def lookup(self, species, stable_id):
    	info = self.perform_rest_action(
    		'/lookup/id/{0}'.format(stable_id),
    		params={'expand':1, 'species': species})
    	return info


@click.group()
@click.pass_context
def main(ctx):
    ctx.obj['client'] = EnsemblRestClient()


@main.command()
@click.argument('seq_id')
@click.option('-t', '--seq_type', default='genomic',
    help='Type of sequence. Defaults to genomic where applicable, i.e. not translations. \
    cdna refers to the spliced transcript sequence with UTR; cds refers to the spliced transcript sequence without UTR.')
@click.pass_context
def get_sequence(ctx, seq_id, seq_type):
	"""
	Get sequence from provided Ensembl ID.
	"""
	sequences = ctx.obj['client'].get_sequence(seq_id, seq_type)
	if sequences:
		print(json.dumps(sequences, sort_keys=True, indent=4))


@main.command()
@click.argument('stable_id')
@click.option('--species', default='human')
@click.pass_context
def lookup(ctx, species, stable_id):
	"""
	Find the species and database for a single identifier e.g. gene, transcript, protein
	"""
	info = ctx.obj['client'].lookup(species, stable_id)
	if info:
		print(json.dumps(info, sort_keys=True, indent=4))

    # client = EnsemblRestClient()
    # variants = client.get_variants(species, symbol)
    # if variants:
    #     for v in variants:
    #         print '{seq_region_name}:{start}-{end}:{strand} ==> {id} ({consequence_type})'.format(**v);

if __name__ == '__main__':
    # if len(sys.argv) == 3:
    #     species, symbol = sys.argv[1:]
    # else:
    #     species, symbol = 'human', 'BRAF'
    # run(species, symbol)
    main(obj={})
