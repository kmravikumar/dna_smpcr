#!/usr/bin/python
import cgi, cgitb, os, sys
cgitb.enable(); # formats errors in HTML

sys.stderr = sys.stdout
print "Content-type: text/html"
print
print """<html>
<head><title>Temperature Conversion</title></head>
<body>
<p>"""

form = cgi.FieldStorage()

if not (form.has_key("from") and form.has_key("value")):
    print "<b>Error</b>: request did not provide the proper query string."
else:
    # This illustrates alternative ways of extracting field
    # parameters. The getfirst method is safer.
    convert = form.getfirst("from")
    v = float(form["value"].value)
    if (convert=='F'):
        C = 5.0*(v-32)/9.0
        print "%.1f F equals %.1f C" % (v,C)
    elif (convert=='C'):
        F = 9.0*v/5.0 + 32
        print "%.1f C equals %.1f F" % (v,F)
    else:
        print "<b>Error</b>: <i>from</i> field must be either F or C"
print "</body></html>"
