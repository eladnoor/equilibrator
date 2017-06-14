eQuilibrator
============

eQuilibrator is a biochemical thermodynamics calculator website.
The website enables calculation of reaction and compound energies
through an intuitive free-text interface. The current online
version of eQuilibrator can be found at:
http://equilibrator.weizmann.ac.il/

If you are looking for the back-end library used to estimate standard Gibbs free energies
[DOI: 10.1371/journal.pcbi.1003098](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003098),
you can find it in the following GitHub repository:
[Component Contribution](https://github.com/eladnoor/component-contribution)

eQuilibrator is written in Python using the Django framework.
It was developed primarily on Ubuntu 16.04 64-bit and is easiest
to develop and set up on that operating system. Setup instructions
below are really only for Ubuntu users.

# Dependencies
- python 2.7
- mysql-server (5.7.16)
- python-glpk
- glpk-tools
- (optional) Indigo Toolkit (https://github.com/ggasoftware/indigo)
Python PyPI:
- Django (1.10.4)
- django-extensions (1.7.5)
- django-haystack (2.5.1)
- glpk (0.4.52)
- matplotlib (1.5.3)
- MySQL-python (1.2.5)
- nltk (3.2.1)
- numpy (1.11.3)
- pandas (0.19.1)
- pip (9.0.1)
- PuLP (1.6.1)
- pyparsing (2.1.10)
- scipy (0.18.1)
- seaborn (0.7.1)
- tablib (0.11.3)
- Solr (3.6.2)
- xlrd (1.0.0)

# Install binary dependencies on Ubuntu
```
sudo apt install mysql-server libmysqlclient-dev python-pip python-dev 
sudo apt install solr-common libglpk-dev glpk-utils
```

# Install missing Python dependencies from PyPI
```
sudo pip install -r requirements.txt
```

# Create MySQL database and user for Django
```
sudo mysql --user=root mysql -p
mysql> CREATE USER '<MYSQLUSR>'@'localhost' IDENTIFIED BY '<MYSQLPWD>';
mysql> GRANT ALL PRIVILEGES ON *.* TO '<MYSQLUSR>'@'localhost';
mysql> CREATE DATABASE <MYSQLDB>;
mysql> exit;
```

# Configure Solr
```
sudo cp solr/schema.xml /etc/solr/conf/
sudo /etc/init.d/tomcat7 restart
```

* Replace the appropriate database name (`<MYSQLDB>`), username (`<MYSQLUSR>`) 
  and password (`<MYSQLPWD>`) in settings.py.
* Run `python manage.py migrate --run-syncdb` to build database schema.
* Run `python db_load_from_sqldump.py` to load the data into the database.
* (optional) instead of `db_load_from_sqldump`, you can use `db_load_from_raw_files`
  which will take much longer (2 hours).

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
sudo cp ~/equilibrator/apache/default.conf /etc/apache2/sites-available/000-default.conf
sudo a2enmod wsgi
sudo apache2ctl restart
```
