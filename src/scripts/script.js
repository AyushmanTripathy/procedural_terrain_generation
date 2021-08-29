const canvas = document.querySelector("canvas");
const ctx = canvas.getContext("2d");

const gridSize = 300;
const blockSize = canvas.width / gridSize;

const pn = new Perlin(Math.random());

const scales = [40, 10];
const persistance = [1, 0.2];
const range = 100;

const frameRate = 20;
let running = false;

start();
function start() {
  if (running) return;
  running = true;
  init(0);
}

function changeScale(index, value) {
  scales[index] = value;
}

function changePersistance(index, value) {
  persistance[index] = value;
}

function stop() {
  running = false;
}

async function init(i) {
  for (let x = i; x < gridSize + i; x++) {
    for (let y = i; y < gridSize + i; y++) {
      //3 layers
      const n1 = Math.round(pn.noise(x / scales[0], y / scales[0], 0) * range);
      const n2 = Math.round(pn.noise(x / scales[1], y / scales[1], 0) * range);

      const n = addLayers(n1, n2);

      ctx.fillStyle = getColor(n);
      ctx.fillRect(
        (x - i) * blockSize,
        (y - i) * blockSize,
        blockSize,
        blockSize
      );
    }
  }

  i += 10;
  await sleep(frameRate);
  console.log("frame");

  if (!running) return;
  ctx.fillStyle = "#202124";
  ctx.fillRect(0, 0, canvas.height, canvas.width);
  init(i);
}

function getColor(num) {
  //snow
  if (num <= 20) return "white";
  //stones
  if (num <= 35) return "#331a00";
  //tree
  if (num <= 40) return "#145214";
  //grass
  if (num <= 60) return "#1f7a1f";
  //sand
  if (num <= 63) return "#eed757";
  //water
  if (num <= 80) return "#0066ff";
  //deep sea
  return "#003380";
}

function addLayers(n1, n2) {
  const sum = n1 * persistance[0] + n2 * persistance[1];
  const max = range * persistance[0] + range * persistance[1];
  return Math.round(map(sum, 0, max, 0, range));
}
function map(n, start1, stop1, start2, stop2) {
  return ((n - start1) / (stop1 - start1)) * (stop2 - start2) + start2;
}

function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}
