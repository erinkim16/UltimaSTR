"""Build a self-contained HTML viewer (str_viewer.html) for the STR result tables.
Open the HTML in any browser: tabs per table, click a header to sort, type in the
boxes under the headers to filter (e.g. 'yes' under lead, 'HTT' under locus,
'>20' or '<0.5' under a number column). Works offline; data is embedded.
Usage: python3 make_viewer.py [table.tsv ...] [--out str_viewer.html]
"""
import argparse, csv, json, os
ap = argparse.ArgumentParser()
ap.add_argument("tables", nargs="*", default=["vincent_str_calls.tiered.tsv", "vincent_str_calls.filtered.tsv", "site_quality.tsv", "cohort_distributions.tsv"])
ap.add_argument("--out", default="str_viewer.html")
a = ap.parse_args()
data = []
for t in a.tables:
    if not os.path.exists(t): print("skip (missing):", t); continue
    rows = list(csv.reader(open(t), delimiter="\t"))
    data.append({"name": os.path.basename(t), "cols": rows[0], "rows": rows[1:]})
html = r"""<!doctype html><html><head><meta charset="utf-8"><title>STR results viewer</title>
<style>
:root{--bg:#fff;--fg:#1d1d1f;--mut:#6b6b70;--line:#e3e3e6;--head:#f5f5f7;--hi:#eaf2ff;--red:#fdeaea;--amb:#fdf3df;--gry:#f1f1f3}
@media (prefers-color-scheme:dark){:root{--bg:#161618;--fg:#ececee;--mut:#9a9aa0;--line:#2d2d31;--head:#202024;--hi:#1e2a3d;--red:#3a1f22;--amb:#3a321e;--gry:#232327}}
body{margin:0;font:13px/1.4 system-ui,-apple-system,Segoe UI,sans-serif;background:var(--bg);color:var(--fg)}
header{padding:12px 16px;border-bottom:1px solid var(--line);display:flex;gap:8px;flex-wrap:wrap;align-items:center}
button.tab{border:1px solid var(--line);background:var(--head);color:var(--fg);padding:6px 12px;border-radius:6px;cursor:pointer}
button.tab.on{background:var(--fg);color:var(--bg)}
#info{color:var(--mut);margin-left:auto}
.wrap{overflow:auto;height:calc(100vh - 58px)}
table{border-collapse:collapse;white-space:nowrap}
th,td{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left}
th{position:sticky;top:0;background:var(--head);cursor:pointer;user-select:none;z-index:2}
tr.f th{top:27px;cursor:default;padding:2px 4px}
tr.f input{width:100%;min-width:60px;box-sizing:border-box;font:inherit;padding:2px 4px;border:1px solid var(--line);border-radius:4px;background:var(--bg);color:var(--fg)}
tbody tr:hover{background:var(--hi)}
tr.path{background:var(--red)} tr.inter{background:var(--amb)} tr.flank td{color:var(--mut)}
.legend span{display:inline-block;padding:1px 6px;border-radius:4px;margin-left:4px}
</style></head><body>
<header><div id="tabs"></div>
<span class="legend"><span style="background:var(--red)">pathogenic range</span><span style="background:var(--amb)">intermediate</span><span style="color:var(--mut)">flanking only</span></span>
<span id="info"></span></header>
<div class="wrap"><table><thead id="th"></thead><tbody id="tb"></tbody></table></div>
<script>
const DATA=__DATA__;let cur=0,sortCol=-1,sortDir=1,filters={};
const num=s=>{const m=String(s).match(/-?\d+(\.\d+)?/);return m?parseFloat(m[0]):NaN};
function tabs(){document.getElementById('tabs').innerHTML=DATA.map((t,i)=>`<button class="tab ${i==cur?'on':''}" onclick="go(${i})">${t.name}</button>`).join(' ')}
function go(i){cur=i;sortCol=-1;filters={};tabs();head();render()}
function head(){const c=DATA[cur].cols;document.getElementById('th').innerHTML=
 '<tr>'+c.map((h,i)=>`<th onclick="srt(${i})">${h}${i==sortCol?(sortDir>0?' ▲':' ▼'):''}</th>`).join('')+'</tr>'+
 '<tr class="f">'+c.map((h,i)=>`<th><input data-i="${i}" placeholder="filter" oninput="filters[${i}]=this.value;render()"></th>`).join('')+'</tr>'}
function srt(i){sortDir=(sortCol==i)?-sortDir:1;sortCol=i;const v=[...document.querySelectorAll('tr.f input')].map(x=>x.value);head();document.querySelectorAll('tr.f input').forEach((x,k)=>x.value=v[k]);render()}
function match(v,f){f=f.trim();if(!f)return true;const m=f.match(/^([<>]=?)\s*(-?[\d.]+)$/);
 if(m){const x=num(v),y=parseFloat(m[2]);if(isNaN(x))return false;return m[1]=='>'?x>y:m[1]=='<'?x<y:m[1]=='>='?x>=y:x<=y}
 return String(v).toLowerCase().includes(f.toLowerCase())}
function render(){const t=DATA[cur],fl=t.cols.indexOf('flag');
 let r=t.rows.filter(row=>Object.entries(filters).every(([i,f])=>match(row[i]??'',f)));
 if(sortCol>=0){r=[...r].sort((a,b)=>{const x=num(a[sortCol]),y=num(b[sortCol]);
  const c=(!isNaN(x)&&!isNaN(y))?x-y:String(a[sortCol]).localeCompare(String(b[sortCol]));return c*sortDir})}
 const show=r.slice(0,2000);
 document.getElementById('tb').innerHTML=show.map(row=>{const f=fl>=0?row[fl]:'';
  const cls=f.includes('PATHOGENIC')?'path':f.includes('INTERMEDIATE')?'inter':f.includes('FLANKING')?'flank':'';
  return `<tr class="${cls}">`+row.map(v=>`<td>${v}</td>`).join('')+'</tr>'}).join('');
 document.getElementById('info').textContent=`${r.length} of ${t.rows.length} rows`+(r.length>2000?' (showing first 2000)':'')}
tabs();head();render();
</script></body></html>"""
open(a.out, "w").write(html.replace("__DATA__", json.dumps(data)))
print("wrote", a.out, "with", ", ".join(f"{d['name']} ({len(d['rows'])} rows)" for d in data))
