eQuilibrator
============

eQuilibrator is a biochemical thermodynamics calculator website. The website enables calculation 
of reaction and compound energies through an intuitive free-text interface. The current online
version of eQuilibrator can be found at:
http://equilibrator.weizmann.ac.il/

eQuilibrator is written in Python using the Django framework. It was developed primarily on Ubuntu
and is easiest to develop and set up on that operating system. Setup instructions below are really
only for Ubuntu users.

# Dependencies
- Django 1.6.2
- MySql 5.5
- Django-Haystack & Xapian (autocomplete)
- PyParsing 1.5.7 (also 2.0.1 tested)
- NLTK 2.0.4
- NumPy and SciPy
- Matplotlib & Seaborn (Plotting)
- Pulp & GLPK (optimization)
- (optional) Indigo Toolkit (https://github.com/ggasoftware/indigo)

# Installing binary dependencies on Ubuntu
```
sudo apt-get install mysql-server libmysqlclient-dev
sudo apt-get install python-pip python-dev
sudo apt-get install python-numpy python-scipy python-matplotlib python-pandas
sudo apt-get install glpk-utils python-glpk
```

Follow this gist to install Xapian 
https://gist.github.com/areski/0919d3b0874fd49ec172

# Remaining Python Dependencies 
```
sudo pip install django==1.7.10 django-extensions django-haystack
sudo pip install seaborn nltk pulp pyparsing MySQL-python
sudo pip install django-debug-toolbar
```

Install a recent version of xapian-haystack from here: 
https://github.com/notanumber/xapian-haystack
```
sudo pip install git+https://github.com/notanumber/xapian-haystack.git
```

# Running the Development Server on a Remote Host

```
sudo python manage.py runserver 0.0.0.0:80
```

Will run the development server on port 80 and allow external IPs to access the server. This is very
useful for debugging differences between your local and remote environments.

# Setting up Apache + Django eQuilibrator
TODO:Finish this.

```
sudo apt-get install apache2 libapache2-mod-wsgi git
git clone https://github.com/eladnoor/equilibrator.git
cd equilibrator
```

* Create a mysql user for Django to use (if needed).
* Create a mysql database for Django to use.
* Put the appropriate database name, username and password in settings.py.
* Edit sqlload.sh to use the correct username and database name. 
* Run ./sqlload.sh to load the database up. 
* Run `python manage.py rebuild_index` to build the xapian index for search.