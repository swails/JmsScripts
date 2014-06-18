#!/usr/bin/python

"""This module provides entry, saveframe, and loop objects. Use python's built in help function for documentation.

There are two variables you can set to control our behavior. Setting bmrb.verbose to True will print some of what is going on to the terminal. Setting raise_parse_warnings to true will raise an exception if the parser encounters something problematic. Normally warnings are suppressed.

Some errors will be detected and exceptions raised, but this does not implement a full validator (at least at present).

Call directly (rather than importing) to run a self-test.
"""

#############################################
#                 Imports                   #
#############################################

import os
import re
import sys
import csv
import bmrb
import copy
import urllib2
import itertools
from cStringIO import StringIO

# Import our clone of an ordered dict if we are using a low version of python
if sys.version_info < (2,7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

from sans import STARLexer
from sans import SansParser
from sans import ErrorHandler, ContentHandler
from sans import quote

#############################################
#            Global Variables               #
#############################################

verbose = False
raise_parse_warnings = False

#############################################
#             Module methods                #
#############################################

def diff(entry1,entry2):
    """Prints the differences between two entries. Non-equal entries will always be detected, but specific differences detected depends on order of entries."""
    diffs = entry1.compare(entry2)
    if len(diffs) == 0:
        print "Identical entries."
    for difference in diffs:
        print difference

def validate(entry,schema=None):
    """Prints a validation report of an entry."""
    validation = entry.validate(schema)
    if len(validation) == 0:
        print "No problems found during validation."
    for err in validation:
        print err

def __cleanValue__(value):
    """Automatically quotes the value in the appropriate way. Don't quote values you send to this method or they will show up in another set of quotes as part of the actual data. E.g.:

    __cleanValue__('"e. coli"') returns '\'"e. coli"\''

    while

    __cleanValue__("e. coli") returns "'e. coli'"

    In fact, you probably will never have to call this method directly. (All print calls automatically use it.)
    """

    if type(value) != str:
        value = str(value)

    # If it's going on it's own line, don't touch it
    if "\n" in value:
        if value[-1] != "\n":
            return value + "\n"
        return value

    if value == "":
        raise ValueError("Empty strings are not allowed as values. Use a '.' or a '?' if needed.")

    # Put tags with newlines on their own ; delimited line
    if " " in value or "\t" in value or value[0] == "_" or "#" in value \
     or (len(value) > 4 and (value[:5] == "data_" or value[:5] == "save_" or value[:5] == "loop_" or value[:5] == "stop_")):
        if '"' in value:
            return "'" + value + "'"
        elif "'" in value:
            return '"' + value + '"'
        elif '"' in value and "'" in value:
            return value + '\n'
        else:
            return "'" + value + "'"

    # It's good to go
    return value

def __formatCategory__(value):
    """Adds a '_' to the front of a tag (if not present) and strips out anything after a '.'"""
    if value:
        if value[:1] != "_":
            value = "_" + value
        if "." in value:
            value = value[:value.index(".")]
    return value

def __formatTag__(value):
    """Strips anything before the '.'"""
    if '.' in value:
        value = value[value.index('.')+1:]
    return value

def __getSchema__(passed_schema=None):
    """If passed a schema (not None) it returns it. If passed none, it checks if the default schema has been initialized. If not initialzed, it initializes it. Then it returns the default schema."""

    global standard_schema
    if passed_schema is None:
        passed_schema = standard_schema
    if passed_schema is None:
        standard_schema = schema()
        passed_schema = standard_schema

    return passed_schema

#############################################
#                Classes                    #
#############################################

class __entryParser__(ContentHandler, ErrorHandler):
    """Parses an entry directly using a SANS parser. You should not ever use this class directly."""
    ent = None
    curframe = None
    curloop = None
    mode = None

    def __init__(self, entry=None):
        if entry is None:
            raise ValueError("You must provide an entry to parse into. Also, why are you using this class?")
        self.ent = entry
        self.curframe = None
        self.curloop = None

    def comment(self, line, text):
        pass

    def startData(self, line, name):
        self.ent.bmrb_id = name

    def endData(self, line, name):
        pass

    def startSaveFrame(self, line, name):
        self.curframe = saveframe.fromScratch(saveframe_name=name)
        self.ent.addSaveframe(self.curframe)

    def endSaveFrame(self, line, name):
        self.curframe = None

    def startLoop(self, line):
        self.curloop = loop.fromScratch()
        self.curframe.addLoop(self.curloop)

    def endLoop(self, line):
        self.curloop = None

    def data(self, tag, tagline, val, valline, delim, inloop):

        if verbose:
            print "Tag / value:", tag, ":", val, "(", tagline, ":", valline, ") d", delim

        if delim == 13:
                val = "$"+str(val)

        if inloop :
            # Update the columns and then add the data
            self.curloop.addColumn(tag,ignore_duplicates=True)
            self.curloop.addDataByColumn(tag,val)
        else :
            self.curframe.addTag(tag,val)

    def error(self, line, msg):
        raise ValueError("Parse error: " + str(msg),line)

    def fatalError(self, line, msg):
        raise ValueError("Fatal parse error: " + str(msg),line)

    def warning(self, line, msg):
        if raise_parse_warnings:
            raise Warning("Parse warning: " + str(msg),line)
        if verbose:
            print "Parse warning: " + str(msg) + " " + str(line)


class schema:
    """A BMRB schema."""

    headers = []
    schema = {}

    def __init__(self, schema_file=None, schema_url=None):
        """Initialize a BMRB schema. With no arguments the most up-to-date schema will be fetched from the BMRB FTP site. Otherwise pass a URL or a file to load a schema from using the schema_url or schema_file optional arguments."""

        self.schema = {}
        self.headers = []
        self.schema_file = None
        self.schema_url = None

        if schema_file is not None and os.path.isfile(schema_file):
            schem_stream = open(schema_file,"rU")
            self.schema_file = schema_file
        else:
            if schema_url is None:
                schema_url = 'ftp://ftp.bmrb.wisc.edu/pub/bmrb/nmr-star_dictionary/dictionary_files/xlschem_ann.csv'
            try:
                with open("/tmp/schema","w") as f:
                    f.write(urllib2.urlopen(schema_url).read())
                schem_stream = open("/tmp/schema","rU")
            except urllib2.HTTPError:
                raise IOError("Could not access STAR schema on BMRB FTP site.")
            self.schema_url = schema_url

        csv_reader = csv.reader(schem_stream)
        self.headers = csv_reader.next()
        for line in csv_reader:
            self.schema[line[8]] = (line[27],line[28],line[1])

    def __repr__(self):
        """Return how we can be initialized."""
        if self.schema_file:
            return "bmrb.schema(schema_file='"+str(self.schema_file)+"')"
        if self.schema_url:
            return "bmrb.schema(schema_url='"+str(self.schema_url)+"')"
        else:
            return "Invalid BMRB schema."

    def __str__(self):
        """Print the schema that we are adhering to."""
        if self.schema_file:
            return "BMRB schema loaded from: '"+str(self.schema_file)+"'"
        elif self.schema_url:
            return "BMRB schema loaded from: '"+str(self.schema_url)+"'"
        else:
            return "BMRB schema from ??? (If you see this you did something wrong.)"

    def valType(self, tag, value, category=None,linenum=None):

        if not tag in self.schema:
            return ["Tag '" + str(tag) + "' not found in schema. Line " + str(linenum) + "."]

        valtype,null_allowed,allowed_category = self.schema[tag]

        if category != None:
            if category != allowed_category:
                return ["The tag '" + str(tag) + "' in category '" + str(category) + "' should be in category '" + str(allowed_category) + "'."]

        if value == ".":
            if null_allowed == "NOT NULL":
                return ["Value cannot be NULL but is: '" + tag + "':'" + value + "' on line " + str(linenum) + "."]
            return []

        if "VARCHAR" in valtype:
            length = valtype[valtype.index("(")+1:valtype.index(")")]
            if len(str(value)) > length:
                return ["Length of value (" + str(len(value)) + ") is too long for VARCHAR(" + length + "): '" + tag + "':'" + value + "' on line " + str(linenum) + "."]
        elif "CHAR" in valtype:
            length = valtype[valtype.index("(")+1:valtype.index(")")]
            if len(str(value)) > length:
                return ["Length of value (" + str(len(value)) + ") is too long for CHAR(" + length + "): '" + tag + "':'" + value + "' on line " + str(linenum) + "."]
        elif "FLOAT" in valtype:
            try:
                a = float(value)
            except Exception as e:
                return ["Value is not of type FLOAT.:'" + tag + "':'" + value + "' on line " + str(linenum) + "."]
        elif "INTEGER" in valtype:
            try:
                a = int(value)
            except Exception as e:
                return ["Value is not of type INTEGER: '" + tag + "':'" + value + "' on line " + str(linenum) + "."]
        return []

class entry:
    """An OO representation of a BMRB entry. You can initialize this object several ways; (e.g. from a file, from the official database, from scratch) see the classmethods."""

    # Put these here for reference
    bmrb_id = 0
    frame_list = []

    def __cmp__(self, other):
        """Returns 1 if this entry is not equal to another entry, 0 if it is equal."""
        return len(self.compare(other))

    def __delitem__(self, item):
        """Remove the indicated saveframe."""

        if isinstance(item,saveframe):
            del self.frame_list[self.frame_list.index(item)]
            return
        else:
            self.__delitem__(self.__getitem__(item))

    def __getitem__(self, item):
        """Get the indicated saveframe."""
        try:
            return self.frame_list[item]
        except TypeError:
            return self.getSaveframeByName(item)

    def __init__(self, **kargs):
        """Don't use this directly, use fromFile, fromScratch, fromString, or fromDatabase to construct."""

        # They initialized us wrong
        if len(kargs) == 0:
            raise ValueError("You must provide either a BMRB ID, a file name, an entry number, or a string to initialize. Use the class methods.")
        elif len(kargs) > 1:
            raise ValueError("You cannot provide multiple optional arguments. Use the class methods instead of initializing directly.")

        # Initialize our local variables
        self.frame_list = []

        if 'the_string' in kargs:
            # Parse from a string by wrapping it in StringIO
            star_buffer = StringIO(kargs['the_string'])
        elif 'file_name' in kargs:
            # Parse directly from a file
            if type(kargs['file_name']) is file:
                star_buffer = kargs['file_name']
            else:
                star_buffer = open(kargs['file_name'], 'r')
        elif 'entry_num' in kargs:
            # Parse from the official BMRB library
            try:
                star_buffer = urllib2.urlopen('http://rest.bmrb.wisc.edu/bmrb/NMR-STAR3/' + str(kargs['entry_num']))
            except urllib2.HTTPError:
                raise IOError("Entry " + str(kargs['entry_num']) + " does not exist in the public database.")
        else:
            # Initialize a blank entry
            self.bmrb_id = kargs['bmrb_id']
            return

        # Load the BMRB entry from the file
        t = __entryParser__(entry=self)
        SansParser.parser( STARLexer(star_buffer), t, t ).parse()

    def __repr__(self):
        """Returns a description of the entry."""
        return "<bmrb.entry '" + str(self.bmrb_id) + "'>"

    def __setitem__(self, key, item):
        """Set the indicated saveframe."""

        # It is a saveframe
        if isinstance(item, saveframe):
            # Add by ordinal
            try:
                self.frame_list[key] = item
            except TypeError:
                # Add by key
                if key in self.frameDict():
                    dict((frame.name,frame) for frame in self.frame_list)
                    for pos,frame in enumerate(self.frame_list):
                        if frame.name == key:
                            self.frame_list[pos] = item
                else:
                    raise KeyError("Saveframe with name '" + str(key) + "' does not exist and therefore cannot be written to. Use the addSaveframe method to add new saveframes.")
        else:
            raise ValueError("You can only assign an entry to a saveframe splice.")

    def __str__(self):
        """Returns the entire entry in STAR format as a string."""
        ret_string = "data_" + str(self.bmrb_id) + "\n\n"
        for frame in self.frame_list:
            ret_string += str(frame) + "\n"
        return ret_string

    @classmethod
    def fromFile(cls, the_file):
        """Create an entry by loading in a file."""
        return cls(file_name=the_file)

    @classmethod
    def fromDatabase(cls, entry_num):
        """Create an entry corresponding to the most up to date entry on the public BMRB server. (Requires ability to initiate outbound HTTP connections.)"""
        return cls(entry_num=entry_num)

    @classmethod
    def fromString(cls, the_string):
        """Create an entry by parsing a string."""
        return cls(the_string=the_string)

    @classmethod
    def fromScratch(cls, bmrb_id):
        """Create an empty entry that you can programatically add to. You must pass a number corresponding to the BMRB ID. If this is not a "real" BMRB entry, use 0 as the BMRB ID."""
        return cls(bmrb_id=bmrb_id)

    def addSaveframe(self, saveframe):
        """Add a saveframe to the entry."""
        self.frame_list.append(saveframe)

    def compare(self, other):
        """Returns the differences between two entries as a list. Otherwise returns 1 if different and 0 if equal. Non-equal entries will always be detected, but specific differences detected depends on order of entries."""
        diffs = []
        if self is other:
            return diffs
        if type(other) is str:
            if str(self) == other:
                return diffs
            else:
                return ['String was not exactly equal to entry.']
        try:
            if str(self.bmrb_id) != str(other.bmrb_id):
                diffs.append("BMRB ID does not match between entries: '"+str(self.bmrb_id)+"' vs '"+str(other.bmrb_id) + "'")
            if len(self.frame_list) != len(other.frame_list):
                diffs.append("The number of saveframes in the entries are not equal: "+str(len(self.frame_list))+" vs "+str(len(other.frame_list)))
            for frame in self.frameDict():
                if other.frameDict().get(frame,None) is None:
                    diffs.append("No saveframe with name '"+str(self.frameDict()[frame].name) + "' in other entry.")
                else:
                    comp = self.frameDict()[frame].compare(other.frameDict()[frame])
                    if len(comp) > 0:
                        diffs.append("Saveframes do not match: '"+str(self.frameDict()[frame].name) + "'")
                        diffs.extend(comp)
        except Exception as e:
            diffs.append("An exception occured while comparing: '" + str(e) + "'")

        return diffs

    def frameDict(self):
        """Returns a dictionary of saveframe name -> saveframe object"""
        return dict((frame.name,frame) for frame in self.frame_list)

    def getLoopsByCategory(self, value):
        """Allows fetching loops by category."""
        results = []
        for frame in self.frame_list:
            for loop in frame.loops:
                if loop.category == __formatCategory__(value):
                    results.append(loop)
        return results

    def getSaveframeByName(self, frame):
        """Allows fetching a saveframe by name."""
        frames = self.frameDict()
        if frame in frames:
            return frames[frame]
        else:
            raise KeyError("No saveframe with name '" + str(frame) + "'")

    def getSaveframesByCategory(self, value):
        """Allows fetching saveframes by category."""
        return self.getSaveframesByTagAndValue("Sf_category", value)

    def getSaveframesByTagAndValue(self, tag_name, value):
        """Allows fetching saveframe(s) by tag and tag value."""

        ret_frames = []

        for frame in self.frame_list:
            try:
                if frame.getTag(tag_name) == value:
                    ret_frames.append(frame)
            except KeyError:
                pass

        return ret_frames

    def printTree(self):
        """Prints a summary, tree style, of the frames and loops in the entry."""
        print repr(self)
        for pos,frame in enumerate(self):
            print "\t[" + str(pos) + "] " + repr(frame)
            for pos2,loop in enumerate(frame):
                print "\t\t[" + str(pos2) + "] " + repr(loop)

    def validate(self,validation_schema=None):
        """Validate an entry against a STAR schema. You can pass your own custom schema if desired, otherwise the schema will be fetched from the BMRB servers. Returns a list of errors found. 0-length list indicates no errors found."""

        errors = []
        for frame in self:
            errors.extend(frame.validate(validation_schema=validation_schema))
        return errors


class saveframe:
    """A saveframe. Use the classmethod fromScratch to create one."""

    tags = OrderedDict()
    loops = []
    name = ""
    tag_prefix = None

    def __cmp__(self, other):
        """Returns 1 if this saveframe is not equal to another saveframe, 0 if it is equal."""
        return len(self.compare(other))

    def __delitem__(self, item):
        """Remove the indicated tag or loop."""

        if isinstance(item,loop):
            del self.loops[self.loops.index(item)]
            return

        else:
            self.__delitem__(self.__getitem__(item))

    def __getitem__(self, item):
        """Get the indicated loop or tag."""
        try:
            return self.loops[item]
        except TypeError:
            try:
                return self.getTag(item)
            except KeyError as e:
                try:
                    return self.loopDict()[item]
                except KeyError:
                    raise KeyError("No tag or loop matching '" + str(item) + "'")

    def __init__(self, **kargs):
        """Don't use this directly. Use the classmethods to construct."""

        # They initialized us wrong
        if len(kargs) == 0:
            raise ValueError("Use the class methods to initialize.")

        # Initialize our local variables
        self.tags = OrderedDict()
        self.loops = []

        if 'the_string' in kargs:
            # Parse from a string by wrapping it in StringIO
            star_buffer = StringIO(kargs['the_string'])
        elif 'file_name' in kargs:
            # Parse directly from a file
            if type(kargs['file_name']) is file:
                star_buffer = kargs['file_name']
            else:
                star_buffer = open(kargs['file_name'], 'r')
        elif 'saveframe_name' in kargs:
            # If they are creating from scratch, just get the saveframe name
            self.name = kargs['saveframe_name']
            if 'tag_prefix' in kargs:
                self.tag_prefix = __formatCategory__(kargs['tag_prefix'])
            return

        star_buffer = StringIO("data_1\n"+star_buffer.read())
        tmp_entry = entry.fromScratch(0)

        # Load the BMRB entry from the file
        t = __entryParser__(entry=tmp_entry)
        SansParser.parser( STARLexer(star_buffer), t, t ).parse()

        # Copy the first parsed saveframe into ourself
        self.tags = tmp_entry[0].tags
        self.loops = tmp_entry[0].loops
        self.name = tmp_entry[0].name
        self.tag_prefix = tmp_entry[0].tag_prefix

    @classmethod
    def fromScratch(cls, saveframe_name, tag_prefix=None):
        """Create an empty saveframe that you can programatically add to. You may also pass the tag prefix as the second argument. If you do not pass the tag prefix it will be set the first time you add a tag."""
        return cls(saveframe_name=saveframe_name, tag_prefix=tag_prefix)

    @classmethod
    def fromFile(cls, the_file):
        """Create a saveframe by loading in a file."""
        return cls(file_name=the_file)

    @classmethod
    def fromString(cls, the_string):
        """Create a saveframe by parsing a string."""
        return cls(the_string=the_string)

    def __repr__(self):
        """Returns a description of the saveframe."""
        return "<bmrb.saveframe '" + str(self.name) + "'>"

    def __setitem__(self, key, item):
        """Set the indicated loop or tag."""

        # It's a loop
        if isinstance(item,loop):
            try:
                integer = int(str(key))
                self.loops[integer] = item
            except ValueError:
                if key in self.loopDict():
                    for pos,tmp_loop in enumerate(self.loops):
                        if tmp_loop.category == key:
                            self.loops[pos] = item
                else:
                    raise KeyError("Loop with category '" + str(key) + "' does not exist and therefore cannot be written to. Use addLoop instead.")
        else:
            # If the tag already exists, set its value
            if key in self.tags:
                self.tags[key] = self.addTag(key,item,ignore_duplicates=True,return_only=True)
            else:
                self.addTag(key,item)

    def __str__(self):
        """Returns the saveframe in STAR format as a string."""

        if self.tag_prefix == None:
            raise ValueError("The tag prefix was never set!")

        # Make sure this isn't a dummy saveframe before proceeding
        try:
            width = max(map(lambda x:len(self.tag_prefix+"."+x), self.tags))
        except ValueError:
            return "\nsave_" + str(self.name) + "\n\nsave_\n"

        # Print the saveframe
        ret_string = "save_" + str(self.name) + "\n"
        pstring = "   %-" + str(width) + "s  %s\n"
        mstring = "   %-" + str(width) + "s\n;\n%s;\n"

        # Print the tags
        for tag in self.tags:
            cleanTag = __cleanValue__(self.tags[tag])
            if "\n" in cleanTag:
                ret_string +=  mstring % (self.tag_prefix+"."+tag, cleanTag)
            else:
                ret_string +=  pstring % (self.tag_prefix+"."+tag, cleanTag)

        # Print any loops
        for loop in self.loops:
            ret_string +=  str(loop)

        # Close 'er up
        ret_string += "save_\n"
        return ret_string

    def addLoop(self, loop):
        """Add a loop to the saveframe loops."""

        if loop.category:
            if loop.category in self.loopDict():
                raise ValueError("You cannot have two loops with the same category in one saveframe. Category: '" + loop.category + "'.")

        self.loops.append(loop)

    def addTag(self, name, value, ignore_duplicates=False, return_only=False):
        """Add a tag to the tag list. Does a bit of validation and parsing. Set ignore_duplicates to true to ignore attempts to add the same tag more than once rather than raise an exception."""

        if "." in name:
            if name[0] != ".":
                prefix = __formatCategory__(name)
                if self.tag_prefix == None:
                    self.tag_prefix = prefix
                elif self.tag_prefix != prefix:
                    raise ValueError("One saveframe cannot have tags with different prefixes (or tags that don't match the set prefix)! " + self.tag_prefix + " vs " + prefix)
                name = name[name.index(".")+1:]
            else:
                name = name[1:]

        # No duplicate tags
        try:
            self.getTag(name)
            if not ignore_duplicates:
                raise ValueError("There is already a tag with the name '" + name + "'.")
            if return_only:
                return value
            else:
                return
        except KeyError:
            pass
        if "." in name:
            raise ValueError("There cannot be more than one '.' in a tag name.")
        if " " in name:
            raise ValueError("Column names can not contain spaces.")

        if verbose:
            print "Adding tag: ("+name+") with value ("+value+")"

        if return_only:
            return value

        self.tags[name] = value

    def compare(self, other):
        """Returns the differences between two saveframes as a list. Non-equal saveframes will always be detected, but specific differences detected depends on order of saveframes."""
        diffs = []

        # Check if this is literally the same object
        if self is other:
            return diffs
        # Check if the other object is our string representation
        if type(other) is str:
            if str(self) == other:
                return diffs
            else:
                return ['String was not exactly equal to saveframe.']
        # Do STAR comparison
        try:
            if str(self.name) != str(other.name):
                diffs.append("\tSaveframe names do not match: "+str(self.name)+"' vs '"+str(other.name) + "'")
            if str(self.tag_prefix) != str(other.tag_prefix):
                diffs.append("\tTag prefix does not match: "+str(self.tag_prefix)+"' vs '"+str(other.tag_prefix) + "'")
            if len(self.tags) != len(other.tags):
                diffs.append("\tNumber of tags does not match: "+str(len(self.tags))+"' vs '"+str(len(other.tags)) + "'")
            for tag in self.tags:
                if self.tags[tag] != other.tags.get(tag,None):
                    diffs.append("\tMismatched tag values for tag '" + str(tag) + "': '" + str(self.tags[tag]).replace("\n","\\n")+"' vs '"+str(other.tags.get(tag,"[not present]")).replace("\n","\\n") + "'")
            if len(self.loops) != len(other.loops):
                diffs.append("\tNumber of children loops does not match: "+str(len(self.loops))+" vs "+str(len(other.loops)))
            for loop in self.loopDict():
                if other.loopDict().get(loop,None) is None:
                    diffs.append("\tNo loop with category '"+str(self.loopDict()[loop].category) + "' in other entry.")
                else:
                    comp = self.loopDict()[loop].compare(other.loopDict()[loop])
                    if len(comp) > 0:
                        diffs.append("\tLoops do not match: '"+str(self.loopDict()[loop].category) + "'")
                        diffs.extend(comp)
        except Exception as e:
            diffs.append("\tAn exception occured while comparing: '" + str(e) + "'")

        return diffs

    def deleteTag(self, tag):
        """Deletes a tag from the saveframe based on tag name."""
        tag = __formatTag__(tag)
        if tag in self.tags:
            del self.tags[tag]
        else:
            raise KeyError("There is no tag with name " + str(tag) + " to remove!")

    def getDataAsCSV(self, header=True, show_category=False):
        """Return the data contained in the loops, properly CSVd, as a string. Set header to False omit the header. Set show_category to True to add the loop category to the headers."""
        csv_buffer = StringIO()
        cwriter = csv.writer(csv_buffer)
        if header:
            if show_category:
                cwriter.writerow([str(self.tag_prefix)+"."+str(x) for x in self.tags])
            else:
                cwriter.writerow([str(x) for x in self.tags])

        data = []
        for piece in [self.tags[x] for x in self.tags]:
            if type(piece) is str and piece[0] == "$":
                piece = piece[1:]
            data.append(piece)

        cwriter.writerow(data)

        csv_buffer.seek(0)
        return csv_buffer.read().replace('\r\n', '\n')

    def getLoopByCategory(self, name):
        """Return a loop based on the loop name (category)."""
        name = __formatCategory__(name)
        for loop in self.loops:
            if loop.category == name:
                return loop
        raise KeyError("No loop with category " + name)

    def getTag(self, query):
        """Allows fetching a tag by name."""

        query = __formatTag__(query)

        if query in self.tags:
            return self.tags[query]
        if self.tag_prefix and (self.tag_prefix + "." + query) in self.tags:
            return self.tags[self.tag_prefix + "." + query]

        raise KeyError("No tag with name '" + str(query) + "'")

    def loopDict(self):
        """Returns a hash of loop category -> loop."""
        return dict((loop.category,loop) for loop in self.loops)

    def loopIterator(self):
        """Returns an iterator for saveframe loops."""
        for loop in self.loops:
            yield loop

    def tagIterator(self):
        """Returns an iterator for saveframe tags."""
        for tag in self.tags:
            yield [tag,self.tags[tag]]

    def printTree(self):
        """Prints a summary, tree style, of the loops in the saveframe."""
        print repr(self)
        for pos,loop in enumerate(self):
            print "\t[" + str(pos) + "] " + repr(loop)

    def validate(self,validation_schema=None):
        """Validate a saveframe against a STAR schema. You can pass your own custom schema if desired, otherwise the schema will be fetched from the BMRB servers. Returns a list of errors found. 0-length list indicates no errors found."""

        # Get the default schema if we are not passed a schema
        my_schema = __getSchema__(validation_schema)

        errors = []
        for pos,tag in enumerate(self.tags):
            errors.extend(my_schema.valType(self.tag_prefix + "." + tag, self.tags[tag], category=self.getTag("Sf_category"), linenum=pos))

        for loop in self.loops:
            errors.extend(loop.validate(validation_schema=validation_schema, category=self.getTag("Sf_category")))

        return errors


class loop:
    """A loop."""

    category = None
    columns = []
    data = []

    def __cmp__(self, other):
        """Returns 1 if this loop is not equal to another loop, 0 if it is equal."""
        return len(self.compare(other))

    def __getitem__(self, item):
        """Get the indicated row from the data array."""
        try:
            return self.data[item]
        except TypeError:
            if type(item) is tuple:
                item = list(item)
            return self.getDataByTag(tags=item)

    def __init__(self, **kargs):
        """Use the classmethods to initialize."""

        # Initialize our local variables
        self.columns = []
        self.data = []
        self.category = None

        # Update our category if provided
        if 'category' in kargs:
            self.category = __formatCategory__(kargs['category'])
            return

        # They initialized us wrong
        if len(kargs) == 0:
            raise ValueError("Use the class methods to initialize.")

        if 'the_string' in kargs:
            # Parse from a string by wrapping it in StringIO
            star_buffer = StringIO(kargs['the_string'])
        elif 'file_name' in kargs:
            # Parse directly from a file
            if type(kargs['file_name']) is file:
                star_buffer = kargs['file_name']
            else:
                star_buffer = open(kargs['file_name'], 'r')
        elif 'saveframe_name' in kargs:
            # If they are creating from scratch, just get the saveframe name
            self.name = kargs['saveframe_name']
            if 'tag_prefix' in kargs:
                self.tag_prefix = __formatCategory__(kargs['tag_prefix'])
            return

        star_buffer = StringIO("data_0\nsave_dummy_frame\n" + star_buffer.read() + "save_")
        tmp_entry = entry.fromScratch(0)

        # Load the BMRB entry from the file
        t = __entryParser__(entry=tmp_entry)
        SansParser.parser( STARLexer(star_buffer), t, t ).parse()

        # Copy the first parsed saveframe into ourself
        self.columns = tmp_entry[0][0].columns
        self.data = tmp_entry[0][0].data
        self.category = tmp_entry[0][0].category


    @classmethod
    def fromScratch(cls, category=None):
        """Create an empty saveframe that you can programatically add to. You may also pass the tag prefix as the second argument. If you do not pass the tag prefix it will be set the first time you add a tag."""
        return cls(category=category)

    @classmethod
    def fromFile(cls, the_file):
        """Create a saveframe by loading in a file."""
        return cls(file_name=the_file)

    @classmethod
    def fromString(cls, the_string):
        """Create a saveframe by parsing a string."""
        return cls(the_string=the_string)


    def __repr__(self):
        """Returns a description of the loop."""
        return "<bmrb.loop '" + str(self.category) + "'>"

    def __str__(self):
        """Returns the loop in STAR format as a string."""

        # If we have no columns than return without printing
        if len(self.columns) == 0:
            if len(self.data) == 0:
                return "\n   loop_\n\n   stop_\n"
            else:
                raise ValueError("Impossible to print data if there are no associated tags. Loop: " + str(self.category))

        # Make sure that if there is data, it is the same width as the column tags
        if len(self.data) > 0:
            for row in self.data:
                if len(self.columns) != len(row):
                    raise ValueError("The number of column tags must match width of the data. Loop: " + str(self.category))

        # Start the loop
        ret_string = "\n   loop_\n"
        # Print the columns
        pstring = "      %-s\n"
        # Check to make sure our category is set
        if self.category == None:
            raise ValueError("The category was never set for this loop. Either add a column with the category intact or specify it when generating the loop.")
        # Print the categories
        for column in self.columns:
            ret_string += pstring % (self.category + "." + column)
        ret_string += "\n"

        if len(self.data) != 0:
            # The nightmare below creates a list of the maximum length of elements in each column in the self.data matrix
            #  Don't try to understand it. It's an imcomprehensible list comprehension.
            title_widths = [max(map(lambda x:len(str(x))+3, col)) for col in [[row[x] for row in self.data] for x in range(0,len(self.data[0]))]]
            # Generate the format string
            pstring = "     " + "%-*s"*len(self.columns) + " \n"

            # Print the data, with the columns sized appropriately
            for datum in copy.deepcopy(self.data):

                # Put quotes as needed on the data
                datum = map(lambda x:__cleanValue__(x), datum)
                for pos,item in enumerate(datum):
                    if "\n" in item:
                        datum[pos] = "\n;\n" + item + ";\n"

                # Print the data (combine the columns widths with their data)
                ret_string += pstring % tuple(itertools.chain.from_iterable([d for d in zip(title_widths,datum)]))

        # Close the loop
        ret_string += "   stop_\n"
        return ret_string

    def addColumn(self, name, ignore_duplicates=False):
        """Add a column to the column list. Does a bit of validation and parsing. Set ignore_duplicates to true to ignore attempts to add the same tag more than once rather than raise an exception.

        You can also pass a list of column names to add more than one column at a time."""

        if isinstance(name, (list, tuple)):
            for x in name:
                self.addColumn(x, ignore_duplicates=ignore_duplicates)
            return

        name = name.strip()

        if "." in name:
            if name[0] != ".":
                category = name[0:name.index(".")]
                if category[:1] != "_":
                    category = "_" + category

                if self.category == None:
                    self.category = category
                elif self.category != category:
                    raise ValueError("One loop cannot have columns with different categories (or columns that don't match the set prefix)!")
                name = name[name.index(".")+1:]
            else:
                name = name[1:]

        # Ignore duplicate tags
        if name in self.columns:
            if ignore_duplicates:
                return
            else:
                raise ValueError("There is already a column with the name '" + name + "'.")
        if "." in name:
            raise ValueError("There cannot be more than one '.' in a tag name.")
        if " " in name:
            raise ValueError("Column names can not contain spaces.")
        self.columns.append(name)

    def addData(self, the_list):
        """Add a list to the data field. Items in list can be any type, they will be converted to string and formatted correctly. The list must be the same cardinality as the column names."""
        if len(the_list) != len(self.columns):
            raise ValueError("The list must have the same number of elements as the number of columns! Insert column names first.")
        # Sanitize the input data and add it
        self.data.append(map(lambda x:str(x), the_list))

    def addDataByColumn(self, column_id, value):
        """Add data to the loop one element at a time, based on column. Useful when adding data from SANS parsers."""
        column_id = __formatTag__(column_id)
        if not column_id in self.columns:
            raise ValueError("The column tag (" + column_id + ") to which you are attempting to add data does not yet exist. Create the columns before adding data.")
        pos = self.columns.index(column_id)
        if len(self.data) == 0:
            self.data.append([])
        if len(self.data[-1]) == len(self.columns):
            self.data.append([])
        if len(self.data[-1]) != pos:
            raise ValueError("You cannot add data out of column order.")
        self.data[-1].append(value)

    def clearData(self):
        """Erases all data in this loop. (But not the data columns.)"""
        self.data = []

    def compare(self, other):
        """Returns the differences between two loops as a list. Otherwise returns 1 if different and 0 if equal. Order of loops being compared does not make a difference on the specific errors detected."""
        diffs = []

        # Check if this is literally the same object
        if self is other:
            return diffs
        # Check if the other object is our string representation
        if type(other) is str:
            if str(self) == other:
                return diffs
            else:
                return ['String was not exactly equal to loop.']
        # Do STAR comparison
        try:
            if str(self.category) != str(other.category):
                diffs.append("\t\tCategory of loops does not match: '"+str(self.category)+"' vs '"+str(other.category) + "'")
            if self.data != other.data or self.columns != other.columns:
                diffs.append("\t\tLoop data does not match for loop with category '"+str(self.category) + "'")
        except Exception as e:
            diffs.append("\t\tAn exception occured while comparing: '" + str(e) + "'")

        return diffs

    def deleteDataByTagValue(self, tag, value, index_tag=None):
        """Deletes all rows which contain the provided value in the provided column. If index_tag is provided, that column is renumbered starting with 1."""

        # Parse the category and tag
        full = str(self.category) + "." + __formatTag__(str(tag))
        tag,category = __formatTag__(full),__formatCategory__(full)

        if category != self.category:
            raise ValueError("Category provided in your column '" + category + "' does not match this loop's category '" + str(self.category) + "'.")

        # Can't delete a row if the tag they want to delete isn't there!
        if tag not in self.columns:
            raise ValueError("The tag you provided '" + tag + "' isn't in this loop!")

        search_column = self.columns.index(tag)

        # Delete all rows in which the user-provided tag matched
        cur_row = 0
        while cur_row < len(self.data):
            if self.data[cur_row][search_column] == value:
                self.data.pop(cur_row)
                continue
            cur_row += 1

        # Re-number if they so desire
        if index_tag is not None:
            self.renumberRows(index_tag)

    def getDataByTag(self, tags=None):
        """Provided a list of tag names, or ordinals corresponding to columns, return the selected tags by row as a list of lists."""

        # All tags
        if tags is None:
            return self.data
        # Turn single elements into lists
        if type(tags) != list:
            tags = [tags]

        # Strip the category if they provide it
        for pos,item in enumerate(map(lambda x:str(x),tags)):
            if "." in item and __formatCategory__(item) != self.category:
                raise ValueError("Cannot fetch data with column (" + str(item) + ") because the category does not match the category of this loop (" + self.category + ").")
            tags[pos] = __formatTag__(item)

        # Map column name to column position in list
        column_mapping = dict(itertools.izip(reversed(self.columns), reversed(xrange(len(self.columns)))))

        # Make sure their fields are actually present in the entry
        column_ids = []
        for query in tags:
            if str(query) in column_mapping:
                column_ids.append(column_mapping[query])
            elif type(query) == int:
                column_ids.append(query)
            else:
                raise ValueError("Your column '" + str(query) + "' was not found in this loop.")

        # Use a list comprehension to pull the correct tags out of the rows
        return [[row[col_id] for col_id in column_ids] for row in self.data]

    def getDataAsCSV(self, header=True, show_category=False):
        """Return the data contained in the loops, properly CSVd, as a string. Set header to False to omit the header. Set show_category to True to add the loop category to the headers."""
        csv_buffer = StringIO()
        cwriter = csv.writer(csv_buffer)

        if header:
            if show_category:
                cwriter.writerow([str(self.category)+"."+str(x) for x in self.columns])
            else:
                cwriter.writerow([str(x) for x in self.columns])

        for row in self.data:

            data = []
            for piece in row:
                if type(piece) is str and piece[0] == "$":
                    piece = piece[1:]
                data.append(piece)

            cwriter.writerow(data)

        csv_buffer.seek(0)
        return csv_buffer.read().replace('\r\n', '\n')

    def printTree(self):
        """Prints a summary, tree style, of the loop."""
        print repr(self)

    def renumberRows(self, index_tag, start_value=1):
        """Renumber a given column incrementally."""

        full = str(self.category) + "." + __formatTag__(str(index_tag))
        tag,category = __formatTag__(full),__formatCategory__(full)

        if tag not in self.columns:
            raise ValueError("The renumbering tag you provided '" + tag + "' isn't in this loop!")

        renumber_column = self.columns.index(tag)

        for x in range(1,len(self.data)+1):
            self.data[x-1][renumber_column] = x

    def validate(self,validation_schema=None, category=None):
        """Validate a loop against a STAR schema. You can pass your own custom schema if desired, otherwise the schema will be fetched from the BMRB servers. Returns a list of errors found. 0-length list indicates no errors found."""

        # Get the default schema if we are not passed a schema
        my_schema = __getSchema__(validation_schema)

        errors = []

        # Check the data
        for rownum,row in enumerate(self.data):
            # Make sure the width matches
            if len(row) != len(self.columns):
                errors.append("Loop '" + str(self.category) + "' data width does not match it's column tag width on row " + str(rownum) + ".")
            for pos,datum in enumerate(row):
                errors.extend(my_schema.valType(self.category + "." + self.columns[pos], datum, category=category, linenum=str(rownum)+" column "+str(pos)))

        return errors




standard_schema = None

# Do unit tests if we are ran directly
if __name__ == '__main__':

    import optparse
    # Specify some basic information about our command
    parser = optparse.OptionParser(usage="usage: %prog",version="%prog",description="STAR handling python module. Usually you'll want to import this. Runs a self test with no arguments.")
    parser.add_option("--diff", action="store_true", dest="diff", default=False, help="Compare two entries.")
    # Options, parse 'em
    (options, cmd_input) = parser.parse_args()

    if options.diff:
        if len(cmd_input) < 2:
            raise ValueError("You must supply two file names as arguments.")
        diff(entry.fromFile(cmd_input[0]), entry.fromFile(cmd_input[1]))
        sys.exit(0)

    print "Running unit tests..."
    new = entry.fromScratch(0)
    new.addSaveframe(saveframe.fromScratch("ham_and_cheese"))
    new.getSaveframeByName('ham_and_cheese').addTag('Sf_category', 'sandwhich')
    new.getSaveframeByName('ham_and_cheese').addTag('_sandwhiches.Flavor', "Too tasty\nfor me\nand you\n")
    new.getSaveframesByCategory('sandwhich')[0].addTag('_sandwhiches.Color', 'Red " and yellow')
    new.getSaveframesByCategory('sandwhich')[0].addTag('Size', '"edible"more or less')
    new.getSaveframesByCategory('sandwhich')[0].addTag('Taste', '"so good"')
    new.getSaveframesByCategory('sandwhich')[0].addTag('Looks', '"so" good')
    new.getSaveframeByName('ham_and_cheese').loops.append(loop.fromScratch())
    t = new.frameDict()['ham_and_cheese']
    t.loops[0].addColumn('_silly.Thing')
    t.loops[0].addData(['goats and sandwhiches'])
    a = new.__str__()

    errors = 0

    def printError(e,x,ent_str):
        if len(e.args) == 1:
            print str(x)+": ",str(e)
        else:
            print str(x)+": ",str(e.args[0]),"on line",(e.args[1]-1)
            splitted = ent_str.split("\n")
            with open("/tmp/"+str(x),"w") as tmp:
                tmp.write(ent_str)
            for x in range(e.args[1]-5,e.args[1]+2):
                try:
                    print "\t%-5d: %s" % (x,splitted[x])
                except IndexError:
                    print "\t%-5d: %s" % (x,"EOF")
                    return

    myrange = (15000,15100)
    if len(sys.argv) == 3: myrange = (int(sys.argv[1]),int(sys.argv[2]))
    if len(sys.argv) == 2: myrange = (int(sys.argv[1]), (int(sys.argv[1])+1))

    use_stardiff = False
    if os.path.exists("/bmrb/linux/bin/stardiff"):
        import subprocess
        use_stardiff = True

    for x in xrange(*myrange):
        try:
            orig_str = urllib2.urlopen('http://rest.bmrb.wisc.edu/bmrb/NMR-STAR3/' + str(x)).read()
        except urllib2.HTTPError:
            continue
        try:
            ent = entry.fromString(orig_str)
            ent_str = str(ent)
        except IOError:
            continue
        except Exception as e:
            printError(e,x,orig_str)
            errors += 1
            continue
        try:
            reent = entry.fromString(ent_str)
            if ent_str != str(reent):
                print str(x)+": Inconsisent output when re-parsed."
                errors += 1
        except Exception as e:
            printError(e,x,ent_str)
            errors += 1
            continue

        if use_stardiff:
            # Write our data to a pipe for stardiff to read from
            with open("/tmp/comparator1","wb") as tmp1:
                tmp1.write(str(ent))
            with open("/tmp/comparator2","wb") as tmp2:
                tmp2.write(str(orig_str))

            compare = subprocess.Popen(["/bmrb/linux/bin/stardiff","-ignore-tag","_Spectral_peak_list.Text_data","/tmp/comparator1","/tmp/comparator2"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)

            # Wait for stardiff to complete
            compare.wait()
            results = compare.stdout.read()
            if not "NO DIFFERENCES REPORTED" in results:
                print str(x)+": Output inconsistent with original: " + results.strip()
                open("/tmp/" + str(x),"wb").write(str(ent_str))
                errors += 1

        comp = ent.compare(reent)
        if len(comp) > 0:
            print str(x)+": Internal entry comparator detects difference(s):"
            diff(ent,reent)

    if errors == 0:
        print "If you didn't see any errors, than everything is working!"
    else:
        print "At least %d errors were found." % (errors)

    sys.exit(0)
