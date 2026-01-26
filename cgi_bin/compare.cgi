#!/usr/bin/env python3
import cgi
import cgitb

cgitb.enable()
print("Content-Type: text/html")
print()

form = cgi.FieldStorage()
query = form.getvalue("query")
comp1 = form.getvalue("comp1")
comp2 = form.getvalue("comp2")

print(query, comp1, comp2)
