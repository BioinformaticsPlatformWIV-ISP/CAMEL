import pandas as pd

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.components.html.htmltableformatter import HtmlTableFormatter
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.scripts.viralconsensuspipeline.snakefile import nextclade3
from camel.scripts.viralconsensuspipeline.snakefile import antivirals

rule antivirals_check_mutations:
    """
    Queries the detected antiviral mutations against the database.
    """
    input:
        DB = '/db/pipelines/viral_consensus/version_1.1/antivirals',
        TSV = lambda wildcards: nextclade3.get_nextclade_output(wildcards, checkpoints, 'TSV', config),
        INFORMS = lambda wildcards: nextclade3.get_informs_subtype(wildcards, checkpoints)
    output:
        JSON = Path(config['working_dir']) / 'antivirals' / 'mutations.json'
    run:
        import json

        # Parse the subtype
        subtype = SnakemakeUtils.load_object(Path(str(input.INFORMS)))['subtype']
        logger.info(f'Subtype: {subtype}')

        # Parse the mutations from the DB
        path_tsv = Path(input.DB, 'mutations.tsv')
        if not path_tsv.exists():
            raise FileNotFoundError(f'Mutations file not found: {path_tsv}')
        data_db_mutations = pd.read_table(path_tsv, keep_default_na=False, na_values='-')
        data_db_mutations = data_db_mutations[data_db_mutations['subtype'] == subtype]
        if len(data_db_mutations) == 0:
            logger.warning(f'No mutations in database for subtype: {len(subtype)}')
        print(data_db_mutations)

        # Parse the associations
        path_tsv_assoc = Path(input.DB, 'associations.tsv')
        if not path_tsv_assoc.exists():
            raise FileNotFoundError(f'Associations file not found: {path_tsv}')
        data_db_associations = pd.read_table(path_tsv_assoc, keep_default_na=False, na_values='-')

        # Parse the detected mutations
        mutations_all = []
        for path_io in [Path(x) for x in input.TSV]:
            path_tsv_nextclade3 = SnakemakeUtils.load_object(path_io)[0].path
            data_n3 = pd.read_table(path_tsv_nextclade3)
            mutations_segment = data_n3.loc[0, 'aaSubstitutions'].split(',')
            print(mutations_segment)
            mutations_segment = [x.split(':')[-1] for x in mutations_segment]
            print(f'Parsed: {len(mutations_segment)} mutations from: {path_tsv_nextclade3}')
            data_segment = pd.DataFrame({
                'mutation': mutations_segment,
                'segment': path_tsv_nextclade3.parent.name.upper()
            })
            mutations_all.extend(data_segment.to_dict('records'))
        data_mutations_detected = pd.DataFrame(mutations_all)
        muts_detected = list(data_mutations_detected[['segment', 'mutation']].itertuples(index=False, name=None))
        print('MUTS DETECTED:')

        def is_present(mut_string: str, d: list[tuple[str, str]], stype: str) -> bool:
            print(mut_string)
            for m in mut_string.split('+'):
                subtype_, segment, mut = m.split('_')
                if subtype_ != stype:
                    return False
                print('  ', segment, mut)
                if mut == 'E119V':
                    import pprint
                    pprint.pprint((segment, mut))
                    import pprint
                    pprint.pprint(d)
                if (segment, mut) not in d:
                    return False
            return True

        # Cross check both
        detected = data_db_mutations.apply(lambda x: (x['segment'], x['mutation']) in muts_detected, axis=1)
        print('MUTS DB:')
        import pprint
        pprint.pprint(list(data_db_mutations[['segment', 'mutation']].itertuples(index=False, name=None)))
        print(f'Detected {sum(detected)} mutations associated with antivirals')
        print(detected)
        print(data_db_mutations[detected])
        keys = set(data_db_mutations[detected]['key'])

        associations = data_db_associations['key'].apply(lambda x: is_present(x, muts_detected, subtype))
        print(f'Detected {sum(associations)} associations with antivirals')
        print(data_db_associations[associations])

        with open(output.JSON, 'w') as handle:
            json.dump({
                'mutations': data_db_mutations[detected].to_dict('records'),
                'associations': data_db_associations[associations].to_dict('records')
            }, handle, indent=2)

rule antivirals_report:
    """
    Creates a HTML report for the antiviral resistance screening.
    """
    input:
        JSON = rules.antivirals_check_mutations.output.JSON
    output:
        HTML = Path(config['working_dir']) / antivirals.OUTPUT_REPORT
    run:
        with open(input.JSON) as handle:
            data_mutations = json.load(handle)

        section = HtmlReportSection('Antiviral resistance')

        # Add table with mutations
        section.add_header('Detected mutations', 4)
        if len(data_mutations['mutations']) == 0:
            section.add_paragraph('No mutations associated with antiviral resistance detected.')
        else:
            columns = [
                {'key': 'subtype', 'title': 'Subtype'},
                {'key': 'segment', 'title': 'Segment'},
                {'key': 'mutation', 'title': 'Mutation'},
            ]
            section.add_table(
                HtmlTableFormatter.format_table_data(pd.DataFrame(data_mutations['mutations']), columns),
                [col['title'] for col in columns], [('class', 'data')]
            )

        # Add table with associations
        section.add_header('Associations', 4)
        if len(data_mutations['associations']) == 0:
            section.add_paragraph('No antiviral resistance associations detected.')
        else:
            def get_color(inh):
                if inh == 'NI':
                    return HtmlTableCell(inh, color='green')
                return HtmlTableCell(inh, color='red')
            print(data_mutations)
            columns = [
                {'key': 'category', 'title': 'Category'},
                {'key': 'key', 'title': 'Key'},
                {'key': 'antiviral', 'title': 'Antiviral'},
                {'key': 'resistance', 'title': 'Resistance', 'fmt': get_color},
            ]
            section.add_table(
                HtmlTableFormatter.format_table_data(pd.DataFrame(data_mutations['associations']), columns),
                [col['title'] for col in columns], [('class', 'data')]
            )

        section.add_header('Additional information', 4)
        section.add_paragraph(
            "The mutations are extracted from the Nextclade output. Abbreviations: normal inhibition (<b>NI</b>), "
            "reduced inhibition (<b>RI</b>), highly reduced inhibition (<b>HRI</b>).")
        SnakemakeUtils.dump_object([ToolIOValue(section)], Path(output.HTML))

rule antivirals_report_empty:
    """
    Creates an empty output report when the assay is disabled.
    """
    output:
        VAL_HTML = Path(config['working_dir']) / antivirals.OUTPUT_REPORT_EMPTY
    run:
        from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
        SnakePipelineUtils.create_empty_report_section('Antiviral resistance', Path(output.VAL_HTML))

rule antivirals_summary:
    """
    Creates the summary output for the antiviral resistance screening.
    """
    input:
        JSON = rules.antivirals_check_mutations.output.JSON
    output:
        TSV = Path(config['working_dir']) / 'antivirals' / 'summary.tsv'
    run:
        with open(input.JSON) as handle:
            data_mutations = json.load(handle)

        data_summary = [
            ('antivirals_mutations', json.dumps(data_mutations['mutations'])),
            ('antivirals_associations', json.dumps(data_mutations['associations'])),
        ]
        with open(output.TSV, 'w') as handle:
            for key, value in data_summary:
                handle.write('\t'.join([key, str(value)]))
                handle.write('\n')
