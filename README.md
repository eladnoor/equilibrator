eQuilibrator
============

eQuilibrator is a biochemical thermodynamics calculator website. The website enables calculation 
of reaction and compound energies through an intuitive free-text interface. The current online
version of eQuilibrator can be found at:
http://equilibrator.weizmann.ac.il/

eQuilibrator is written in Python using the Django framework. It was developed primarily on Ubuntu 16.04 64-bit and is easiest to develop and set up on that operating system. Setup instructions below are really only for Ubuntu users.

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
- Whoosh (2.7.4)
- xlrd (1.0.0)

# Install binary dependencies on Ubuntu
```
sudo apt install mysql-server libmysqlclient-dev python-pip python-dev glpk-tools python-glpk
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

* Replace the appropriate database name (<MYSQLDB>), username (<MYSQLUSR>) 
  and password (<MYSQLPWD>) in settings.py.
* If necessary, run `python manage.py reset_db` to clear history.
* Run `python sql_load.py` to load the database up. 
* Run `python manage.py rebuild_index` to build the search index.

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
