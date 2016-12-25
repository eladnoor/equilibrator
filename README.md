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
- Django 1.7.10
- Django extensions 1.5.7
- Django profiler 2.0.0
- Django haystack 2.4.0
- Xapian 1.2.23 (indexing for search)
- Xapain haystack 2.0.0
- MySql 5.7.16
- PyParsing 2.1.8
- NLTK 2.0.4
- NumPy 1.11.1
- SciPy 0.18.1
- Matplotlib & Seaborn (Plotting)
- Pulp & GLPK (optimization)
- (optional) Indigo Toolkit (https://github.com/ggasoftware/indigo)

# Install binary dependencies on Ubuntu
```
sudo apt install mysql-server libmysqlclient-dev
sudo apt install python-pip python-dev
sudo apt install python-numpy python-scipy python-matplotlib python-pandas
sudo apt install glpk-utils python-glpk uuid-dev zlib1g-dev g++
sudo apt install python-xapian libxapian-dev
```

# Install Python dependencies from PyPI
```
sudo pip install django==1.7.10
sudo pip install django-extensions==1.5.7
sudo pip install xapian-haystack==2.0.0
sudo pip install django-haystack==2.4.0
sudo pip install seaborn nltk pulp pyparsing MySQL-python tablib
```

# Create MySQL database and user for Django
```
sudo mysql --user=root mysql -p
mysql> CREATE USER 'milolab_eqbtr'@'localhost' IDENTIFIED BY '******';
mysql> GRANT ALL PRIVILEGES ON *.* TO 'milolab_eqbtr'@'localhost';
mysql> CREATE DATABASE milolab_eqbtr;
mysql> exit;
```

* Put the appropriate database name, username and password in settings.py.
* Edit sqlload.sh to use the correct username and database name. 
* Run `./sqlload.sh` to load the database up. 
* Run `python manage.py rebuild_index` to build the xapian index for search.

# Running the Development Server on a Remote Host

```
sudo python manage.py runserver 0.0.0.0:80
```

Will run the development server on port 80 and allow external IPs to access the server. This is very
useful for debugging differences between your local and remote environments.

# Setting up Apache + Django eQuilibrator

* For running eQuilibrator on a web server that can be accessed safely from the internet.
```
sudo apt install apache2 libapache2-mod-wsgi links
git clone https://github.com/eladnoor/equilibrator.git
cd equilibrator
chmod 666 gibbs.log
sudo cp ~/equilibrator/apache2/default.conf /etc/apache2/sites-available/000-default.conf
sudo a2enmod wsgi
sudo apache2ctl restart
```

