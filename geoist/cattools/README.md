# CatalogueTool-Lite
*A simplified Python toolkit for earthquake catalogue manipulation*

The toolkit consists of 9 main modules:
  * **Catalogue** - The main module for database manipulation an I/O
  * **Parsers** - Ad-hoc parsers for specific catalogues and bulletins (ISC, GCMT...)
  * **Selection** - Functions for high-level manipulation of catalogue objects
  * **Exploration** - Functions to explore database information and perform basic statistical analysis
  * **Regressor** - Utilities for magnitude conversion and catalogue homogenisation
  * **MagRules** - Library of magnitude conversion functions
  * **MapTools** - Utility to plot earthquake databases on a map

<!--- -------------------------------------------------------------------------------- --> 
## 1 - The Module *Catalogue*

### 1.1 - Database Initialisation
A catalogue object can be instantiated using the *Database* class from the *Catalogue* module as:
~~~
import Catalogue as Cat
Db0 = Cat.Database('Some Name (Optional)','Some Information (Optional)')
~~~
Optional parameters are the catalogue name and a description string.

### 1.2 - Object Structure

#### 1.2.1 - Attributes
The most important attributes of the catalogue objects are the *Header* and *Events* variables. While *Header* is basically just a dictionary storing diverse information about the catalogue (e.g. name, some descrition...), *Events* is the actual database (as list) of earthquake records, with a more complex internal structure.
Each element of the *Events* list is itself a dictionary containing data of a single event, grouped in four main keys: *Id*, *Magnitude*, *Location* and *Log*.
Here it is an example for one event:
~~~python
Db.Events[0] = {'Id': '02000',
                'Location': [{'LocCode': 'ISC',  # First Location Solution
                              'Year': 1972,
                              'Month': 1,
                              'Day': 19,
                              'Hour': 0,
                              'Minute': 37,
                              'Second': 7.5,
                              'Longitude': -13.8138,
                              'Latitude': 31.3611,
                              'Depth': 33.0,
                              'LatError': None,
                              'LonError': None,
                              'DepError': None,
                              'SecError': 0.35,
                              'Prime': True}],
                'Magnitude': [{'MagCode': 'EMEC',  # First Magnitude Solution
                               'MagSize': 5.0,
                               'MagError': None,
                               'MagType': 'Mw'},
                               {'MagCode': 'ISC',  # Second Magnitude Solution
                               'MagSize': 4.9,
                               'MagError': 0.1,
                               'MagType': 'Ms'}],
                'Log': 'MERGED(EMEC Africa:1234);PREID(1111);'}
~~~
*Location* and *Magnitude* are also list of dictionaries. Each elements of those lists represents a separate solution from a particular agency. In the example above, the event *02000* contains two independent magnitude solutions, but only one location solution.
*Log* is jut a container for processing information (explained later), although it could be used to store any arbitrary text data of an event.

#### 1.2.2 - Methods
Several methods for database manipulation, I/O and exploration are available:
  * *AddEvent* - Add an earthquake event to the database
  * *DelEvent* - Remove an earthquake avent from the database
  * *Import* - Import catalogue from file (csv format)
  * *Export* - Export catalogue to file (csv format)
  * *Load* - Import database structure from binary file (cPickle compressed)
  * *Dump* - Exprot database structure to binary file (cPickle compressed)
  * *Filter* - Filter earthquake events by key field and rule
  * *SetField* - Set database key field to a specific value
  * *Extract* - Extract database information by key field
  * *KeyStat* - Compute statistics on key field occurrence
  * *Print* - Print event information on screen (by ID or index)
  * *Copy* - Create hardcopy of the database
  * *Append* - Concatenate event list of two databases
  * *Size* - Output number of earthquake events
  * *Sort* - Sort events according to origin time
  * *GetIndex* - Get event index from ID string
  * *SetID* - Regenerate progressive IDs

In the following we describe the the use of the different methods, grouped by category.

### 1.3 - Catalogue I/O
Once instantiated an catalogue object, the database can be inflated manually (element by element) or by parsing an external source file. A parsed catalogue, however, can also be manually augmented with new available information.

#### 1.3.1 - Creating a Catalogue Manually
As an example, database items can be created manually in this way:
~~~python
import Catalogue as Cat
Db = Cat.Database('MyCat')

L = [{'Year': 1961, 'Month': 12, 'Day': 3, 'Hour': 5, 'Minute': 20, 'Second': 10}},
     {'Year': 1962, 'Month': 2, 'Day': 5, 'Hour': 12, 'Minute': 6, 'Second': 5}]

M = [{'MagCode': 'AAA', 'MagSize':5, 'MagError': 0.1, 'MagType':'Mw'},
     {'MagCode': 'BBB', 'MagSize':7, 'MagError': 0.2, 'MagType':'ML'}]

# Creating an new empty catalogue item
Db.AddEvent('E001')

# Creating a new item with just Location information
Db.AddEvent('E002', Location=L)

# Creating a new item with Location and Magnitude information
Db.AddEvent('E003', Location=L, Magnitude=M)

# Adding new information to an existing item
Db.AddEvent('E001', L, [], Append=True)
Db.AddEvent('E002', Magnitude=M, Append=True)

# Remove an existing item (by ID)
Db.DelEvent('E003')
~~~
Providing all key fields is however not compulsory. The *AddEvent* method will include only the available information.

#### 1.3.2 - Reading/Writing ASCII Files
Manual creation of a large number of items is however impractical. To avoid that, an earthquake catalogue can be parsed from a csv (ascii) file. The standard format used in the toolkit is:

| Id | Year | Month | Day | Hour | Minute | Second | Longitude | Latitude | Depth | DepError | LocCode | MagSize | MagError | MagType | MagCode |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 16957871 | 1905 | 9 | 8 | 1 | 43 | 2.24 | 15.7840 | 38.6360 | 15.00 | 6.70 | ISC-GEM | 7.16 | 0.70 | Mw | ISC-GEM |
| 16957879 | 1905 | 12 | 4 | 12 | 20 | 7.96 | 38.8300 | 37.2160 | 15.00 | 5.60 | ISC-GEM | 5.56 | 0.63 | Mw | ISC-GEM |
| 16958009 | 1908 | 12 | 28 | 4 | 20 | 26.62 | 15.3490 | 38.0000 | 15.00 | 6.00 | ISC-GEM | 7.03 | 0.37 | Mw | ISC-GEM |

Keywords are the same as those used in the *Events* structure.

A standard csv catalogue can then be parsed using the method *Import*:
~~~python
Db.Import('data/isc-rev-africa-select.csv')
~~~

Nonetheless, non-standard csv formats can also be parsed using *Import*. In this case, the file format (either column names and fields to skip) must be explicitly specified through the *Header* variable. For example:
~~~python
H = ['Id','','Year','Month','Day','Hour','Minute','Second',
     'Longitude','Latitude','','','','Depth','DepError',
     'MagSize','MagError','','','','','','','','','']

Db = Cat.Database('ISC-GEM')
Db.Import('data/isc-gem-v3.csv', Header=H, SkipLine=1, Delimiter=',')
~~~

Identically, the catalogue object can be exported in csv standard format with:
~~~python
Db.Export('data/isc-rev-africa-select.csv')
~~~
Unfortunately, a significant limitation of using the standard CATK format is that only one solution (magnitude or location) is possible per event, either when importing or exporting a csv file. To circumvent the problem, ISF format (see *Parser* module) and binary I/O (next section) should be used instead.

#### 1.3.3 - Reading/Writing Binary Files
To speedup disk I/O access of large catalogue objects, binary (CPickle) files can be used instead of csv. This can simply be done with:
~~~python
# Writing to binary file
Db.Dump('data/isc-gem-v3.bin')
# Reading from binary file
Db.Load('data/isc-gem-v3.bin')
~~~

### 1.4 - Basic Catalogue Manipulation

#### 1.4.1 - Event Selection by Key/Value
Probably, the most useful method for event selection is *Filter*, which allows removing events from a catalogue according to user-defined rules. The method operates on a given Magnitude or Location key by filtering out all non-matching events.

By default, an equality match is performed:
~~~python
# Keep only events with one or more ISC Location solutions
Db.Filter('LocCode', 'ISC')
# Multiple values can also be provided
Db.Filter('LocCode', ['ISC','GCMT','NEIC'])
~~~
However, inequality matching is also allowed, by specifying the argument *Opr*:
~~~python
# Select magnitude equal or higher than 5
Db.Filter('MagSize', 5, Opr='>=')
# Remove unknown depth solutions
Db.Filter('Depth', None, Opr='!=')
~~~
The optional boolean arguments *All* and *Best* can be used to perform conditional selection of multiple entries. *All* is used to select items which contain all values in the match list (equivalent to boolean `and`). This is often the case when comparing multiple solutions from different agencies. *Best*, on the contrary, selects items according to the occurence of just the first matching element of the list (equivalent to boolean `or`). This is very useful for agency prioritisation.
~~~python
Db.Filter('LocCode', ['ISC','GCMT'], All=True)
Db.Filter('LocCode', ['ISC','GCMT'], Best=True)
~~~
By default, the *Filter* method remove non-matching items from the database object permanently. To avoid this behaviour, an hard-copy of the object can be created in output by setting the boolean argument *Owrite* to `False`:
~~~python
DbNew = Db.Filter('LocCode', ['ISC','GCMT'], Owrite=True)
~~~

#### 1.4.2 - Modifying Database Fields
Each database field can be individually accessed/modified through the *Events* structure. The method *SetField* is used to modify simultaneously all entries of the database for a given Key.
~~~python
Db.SetField('MagCode', 'GEM')
~~~
The optional argument *Match* can be used to modify only those entries matching a specific key/value pair (conditional filtering).
~~~python
Db.SetField('MagType', 'Mw', Match=['MagCode','GEM'])
~~~

#### 1.4.3 - Sorting Events Chronologically
After merging information from separate databases, event list might not be sorted chronologically anymore. The method *Sort* can be used to sort items within a second precision.
~~~python
Db.Sort()
~~~
In combination, the method *SetID* can be used to regenerate events ID label according to a progressive numbering. Optional formatting strings can be added before and/or after the ID number:
~~~python
Db.SetID(Str0='Before_number', Str1='After_Number'):
~~~

#### 1.4.4 - Copying and Merging Databases
A hard-copy of a whole database object can be performed using the method *Copy*:
~~~python
DbNew = Db.Copy()
~~~
Events from one catalogue can be appended to the event list of a second one by using the method *Append*:
~~~python
Db.Append(Db2)
~~~
This method, however does not search for or merges duplicated events (for that, we refer to the module *Selection*).

### 1.5 - Extracting Database Information
The method *Size* provides the number of items (events) contained in a database object. Single event information can then be obtained by list index of ID string using the method *Print*:
~~~python
Db0.Print(2000)
~~~
with output:
~~~shell
Event Id: 02000
Location:
[0] - Year: 1972 Month: 2 Day: 1 Hour: 11 Minute: 42 Second: 23.15 Latitude: 35.3609 Longitude: -4.5653 Depth: 41.6 Agency: ISC Prime: True
Magnitude:
[0] - Type: Mw Size: 4.09 Error: 0.28 Agency: NEIS
Log:
MAGCONV(NEIS:mb);PREID(775874);PREID(02004);
~~~

<!--- -------------------------------------------------------------------------------- --> 
## 2 - The Module *Parsers*

The module *Parsers* extends the functionalities of the standard database object to import fixed-format catalogues and bulletins.
Presently available parsers are for the ISC bulletin in ISF format:
~~~python
import Parsers as Par
Db = Par.Database('ISC-AFRICA')
Db.ImportIsf('data/isc-rev-africa.isf')
~~~
and for GCMT catalogue in NDK format:
~~~python
import Parsers as Par
Db = Par.Database('GCMT')
Db.ImportNdk('data/jan76_dec13.ndk')
~~~
