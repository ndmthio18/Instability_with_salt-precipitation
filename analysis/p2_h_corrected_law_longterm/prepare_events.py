from pathlib import Path
import csv
R=Path(__file__).resolve().parent
with (R/'cases_p1.csv').open() as f:
    reader=csv.DictReader(f);fields=reader.fieldnames
    rows=[row for row in reader if 'beta05' in row['id'] or 'av16' in row['id']]
for row in rows:row['id']+='_event'
with (R/'cases_events.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fields,lineterminator='\n');writer.writeheader();writer.writerows(rows)
