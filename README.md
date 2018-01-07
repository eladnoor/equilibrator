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
[Component Contribution](https://github.com/eladnoor/equilibrator-api)

eQuilibrator is written in Python using the Django framework.
It was developed primarily on Ubuntu 16.04 64-bit and is easiest
to develop and set up on that operating system. Setup instructions
below are really only for Ubuntu users.

# Ubuntu dependencies (apt):
- python 3.5+
- mysql-server (5.7.20)
- libmysqlclient-dev (5.7.20)
- glpk-utils (4.63)
- libglpk-dev (4.63)
- solr-tomcat (3.6.2)
- (optional) Indigo Toolkit (https://github.com/ggasoftware/indigo)
# Python PyPI dependencies (pip):
- Django (2.0)
- django-extensions (1.9.7)
- django-haystack (2.6.1)
- matplotlib (2.1.0)
- mysqlclient (1.3.12)
- nltk (3.2.5)
- numpy (1.13.3)
- pandas (0.21.0)
- pulp (1.6.8)
- pyparsing (2.2.0)
- scipy (1.0.0)
- seaborn (0.8.1)
- tablib (0.12.1)
- solr (0.4)
- xlrd (1.1.0)

# Install Ubuntu dependencies
```
sudo apt install mysql-server libmysqlclient-dev python3-pip python3-dev solr-tomcat libglpk-dev glpk-utils
```

# Install Python dependencies
```
sudo pip install -r requirements.txt
```

# Configure Solr
```
sudo cp solr/schema.xml /etc/solr/conf/
sudo /etc/init.d/tomcat8 restart
```
* Make sure solr is running by going to http://127.0.0.1:8080/solr/

# Create MySQL database and synchronize the data
```
python init_db.py
```

# Running the Development Server on a Remote Host

```
sudo python manage.py runserver 0.0.0.0:8000
```

Will run the development server on port 8000 and allow external IPs to access the server. This is very
useful for debugging differences between your local and remote environments.

# Setting up Apache + Django eQuilibrator

For running eQuilibrator on a web server that can be accessed safely from the internet.
Open the file equilibrator/settings.py for editing:
* add the IP of the server to the ALLOWED_HOSTS
* change DEBUG to False
* run the following from the command line:
```
sudo apt install apache2 libapache2-mod-wsgi-py3 links
git clone https://github.com/eladnoor/equilibrator.git
cd equilibrator
sudo python3 manage.py collectstatic
sudo cp apache/default.conf /etc/apache2/sites-available/000-default.conf
sudo a2enmod wsgi
sudo apache2ctl restart
```

