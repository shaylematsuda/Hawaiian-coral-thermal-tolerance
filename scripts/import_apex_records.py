"""Import supplied raw APEX files without changing them.
Exact duplicates are collapsed, differing values at a tank/timestamp remain
unresolved and are written to an audit table rather than averaged.
The headerless 14marl17 tab inherits tanks 1-4 from its source workbook.
Other tabs use their actual Tmp.* headers, including the June workbook's 5-8 tab.
"""
from pathlib import Path
import openpyxl,csv,datetime as dt,collections,math
ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'Data/APEX/derived'
out.mkdir(parents=True,exist_ok=True)
data=collections.defaultdict(list)
def parse(x):
 if isinstance(x,dt.datetime):return x
 if not isinstance(x,str):return None
 for fmt in ('%H:%M %m/%d/%Y','%m/%d/%y %H:%M','%m/%d/%Y %H:%M','%Y-%m-%d %H:%M:%S'):
  try:return dt.datetime.strptime(x.strip(),fmt)
  except ValueError:pass
 return None
files=['14march2017 5-8.xlsx','14March2017 Tank9-12.xlsx','APEX Temp Tanks 1-4.xlsx','Tanks 1-4 June2017.xlsx']
for name in files:
 w=openpyxl.load_workbook(ROOT/'Data/APEX/source_workbooks'/name,read_only=True,data_only=True)
 for s in w:
  headers=None
  for rownum,row in enumerate(s.iter_rows(values_only=True),1):
   if not row or not any(x is not None for x in row):continue
   if row[0]=='Date':headers=list(row);continue
   date=parse(row[0])
   if date is None:continue
   if headers is None:
    if name=='APEX Temp Tanks 1-4.xlsx':headers=['Date','Tmp.1','Tmp.2','Tmp.3','Tmp.4']
    else:continue
   for i,key in enumerate(headers):
    if isinstance(key,str) and key.startswith('Tmp.') and i<len(row) and isinstance(row[i],(int,float)) and math.isfinite(row[i]):
     data[(int(key.split('.')[1]),date)].append((round(float(row[i]),6),name,s.title,rownum))
for p in (ROOT/'Data/APEX').glob('Tank*.csv'):
 with p.open(newline='') as f:
  for rn,row in enumerate(csv.DictReader(f),2):
   date=parse(row['Date'])
   for k,v in row.items():
    if k.startswith('Tmp.'):
     try:value=float(v)
     except (ValueError,TypeError):continue
     data[(int(k.split('.')[1]),date)].append((round(value,6),p.name,'CSV',rn))
records=[];conflicts=[];details=[];gaps=[];summary=[]
for tank in sorted({k[0] for k in data}):
 segment=0;prev=None;nc=0;maxspread=0
 for key in sorted((k for k in data if k[0]==tank),key=lambda k:k[1]):
  date=key[1]; values=sorted({x[0] for x in data[key]});isconflict=len(values)>1
  if prev and (date-prev).total_seconds()>3600:
   segment+=1
   if (date-prev).total_seconds()>86400:gaps.append([tank,prev.isoformat(' '),date.isoformat(' '),(date-prev).total_seconds()/3600])
  if isconflict:
   segment+=1;nc+=1;spread=values[-1]-values[0];maxspread=max(maxspread,spread)
   conflicts.append([tank,date.isoformat(' '),values[0],values[-1],spread,len(values)])
   for value,file,sheet,row in data[key]:details.append([tank,date.isoformat(' '),value,file,sheet,row])
  else:records.append([tank,date.isoformat(' '),values[0],segment])
  if isconflict:segment+=1
  prev=date
 summary.append([tank,nc,maxspread])
def write(name,headers,rows):
 with (out/name).open('w',newline='') as f:w=csv.writer(f);w.writerow(headers);w.writerows(rows)
write('agreeing_records.csv',['Tank','Date','Temperature','Segment'],records)
write('conflicting_timestamps.csv',['Tank','Date','Min','Max','Difference_C','Distinct_values'],conflicts)
write('conflict_source_rows.csv',['Tank','Date','Temperature','File','Sheet','Row'],details)
write('gaps_over_24h.csv',['Tank','Start','End','Hours'],gaps)
write('summary.csv',['Tank','Conflicting_timestamps','Max_difference_C'],summary)
print('Conflicts by tank:',summary)
print('Conflict dates:',min(x[1] for x in conflicts),max(x[1] for x in conflicts))
print('Largest differences:',sorted(conflicts,key=lambda x:-x[4])[:5])
print('Internal gaps >24h:',gaps)
