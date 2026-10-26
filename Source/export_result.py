def save_result(input_path, output_path, original_record, result_record):
    from pathlib import Path
    from io import StringIO
    from Bio import SeqIO
    from dna_writer import write_annotated_dna
    
    input_p = Path(input_path)
    output_p = Path(output_path)
    
    if output_p.suffix.lower() == '.dna':
        if not input_p.suffix.lower() == '.dna':
            raise ValueError('Input must be .dna when output is .dna')
        if result_record.seq != original_record.seq:
            raise ValueError('Sequence mismatch for .dna output')
        if result_record.features[:len(original_record.features)] != original_record.features:
            raise ValueError('Features mismatch for .dna output')
        write_annotated_dna(input_p, output_p, result_record.features[len(original_record.features):])
    elif output_p.suffix.lower() in ('.gb', '.gbk', '.genbank'):
        buffer = StringIO()
        SeqIO.write(result_record, buffer, 'genbank')
        with open(output_p, 'x', encoding='utf-8') as f:
            f.write(buffer.getvalue())
    else:
        raise ValueError(f'Unsupported output suffix: {output_p.suffix}')
