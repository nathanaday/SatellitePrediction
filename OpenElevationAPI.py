import requests

lat = 41.161758
lon = -8.583933
url = r"https://api.open-elevation.com/api/v1/lookup?locations=41.161758,-8.583933"

res = requests.get(url)
print(res.json())
