#!/usr/bin/env python3
import cgi
import cgitb

cgitb.enable()
print("Content-Type: text/html")
print()
print(
    '<form action = "/cgi-bin/compare.cgi" method = get> Enter cancer project (e.g. TCGA-LUAD): <input type = "text" name = "query"> <br />Compound1: <input type = text name = comp1 /> <br />Compound2: <input type = text name = comp2 /> <input type = submit value = Submit /> </form>'
)
