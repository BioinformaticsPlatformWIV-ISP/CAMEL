from app.components.html.htmlreport import HtmlReport


rule initialize_report:
    """
    Initializes the HTML report.
    Creates a HTML file with the title, input section and parameter section.
    """
    output:
        HTML=config.get('report')
    params:
        assembler = config['assembler'],
        fastq_pe = config['fastq_pe'],
        sample_name = config['sample_name'],
        detection_method = config['detection_method']
    run:
        from resources import CSS_STYLE
        import datetime
        report = HtmlReport(output.HTML)
        report.initialize('VTEC Pipeline', CSS_STYLE)
        report.add_header('VTEC Pipeline on {}'.format(params.sample_name), 1, [('class', 'top_header')])

        report.add_header('Input', 2)
        report.add_table([
            ['Input files:', ', '.join(os.path.basename(f) for f in params.fastq_pe)],
            ['Analysis date:', datetime.datetime.now().strftime('%d/%m/%Y - %X')],
            ['Pipeline version:', __PIPELINE_VERSION]],
            table_attributes=[('class', 'information')])

        report.add_header('Settings', 2)
        report.add_table([
            ['Assembler:', params.assembler],
            ['Detection method:', params.detection_method]],
            table_attributes=[('class', 'information')])
        report.save()
