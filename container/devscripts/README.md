# devscripts
development scripts for image maintenance

# notes

### building docker images
- user must be in `docker` group, `pwd` must contain `docker-compose.yml`
```
docker-compose build dev
```

### pushing to dockerhub
example for CI image: give it the dockerhub repo name, with the tag `dev`
```
docker login
docker images
docker tag largexeic_dev:latest cjdilks/largex-eic:dev
docker scan cjdilks/largex-eic:dev  2>&1 | tee vulnerabilities.txt
docker push cjdilks/largex-eic:dev
docker logout
```
