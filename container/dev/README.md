# devscripts
development scripts for image maintenance

# notes

### building docker images
- user must be in `docker` group, `pwd` must contain `docker-compose.yml`
```
docker-compose build dev
container/docker-shell.sh   # for testing
```

### pushing to dockerhub
```
docker login
docker images
docker tag sidiseic_dev:latest cjdilks/sidis-eic:latest
docker scan cjdilks/sidis-eic:latest  2>&1 | tee vulnerabilities.txt
docker push cjdilks/sidis-eic:latest
docker logout
# clean ~/.docker/config.json (if necessary)
```
