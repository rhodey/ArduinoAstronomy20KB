const mini       = require('minimist')
const SerialPort = require('serialport')
const Readline   = require('@serialport/parser-readline')
const model      = require('./model.js')

function onClose(err) {
  if (err) {
    console.error(err)
    process.exit(1)
  } else {
    process.exit(0)
  }
}

function closeIfErr(err) {
  if (err) onClose(err)
}

function sleep(ms, val) {
  return new Promise(function (res, rej) {
    setTimeout(() => res(val), ms)
  })
}

function openSerial(path) {
  return new Promise(function (res, rej) {
    const serial = new SerialPort(path, { baudRate: 115200 })
    serial.once('open', () => res(serial))
    serial.once('error', rej)
  })
}

function Writer(serial) {
  const parser = new Readline()
  serial.pipe(parser)
  let callback = null
  let retry = true

  function read(line) {
    if (line.indexOf('unexpected') >= 0) {
      onClose(line)
    } else if (line.indexOf('error') < 0) {
      retry = true
      callback()
    } else if (retry) {
      retry = false
      console.log(line)
      callback(true)
    } else {
      onClose(line)
    }
  }

  parser.on('data', read)

  return function w(write) {
    return new Promise(function (res, rej) {
      let timer = setTimeout(() => rej(new Error(`write ${write.address} timed out`)), 2000)
      callback = function (retry) {
        clearTimeout(timer)
        if (!retry) {
          res(write.address)
        } else {
          console.log(`retry address ${write.address} ${write.type} ${write.val}.`)
          w(write)
        }
      }
      serial.write(`${write.address}${write.type}${write.val}\n`, function (err) {
        if (err) { rej(err) }
      })
    })
  }
}

const argv = mini(process.argv.slice(2))
if (argv._.length <= 0) { onClose(new Error('node lib/js/flash.js <path>')) }
process.on('SIGINT', onClose)

openSerial(argv._[0])
  .then(serial => sleep(2000, serial))
  .then(function (serial) {
    console.log('serial opened.')
    const write = Writer(serial)
    const writes = model()
    let i = 0

    function next() {
      if (i >= writes.length) {
        console.log(`end address ${writes[writes.length - 1].address + 4}.`)
        return serial.close()
      } else {
        return write(writes[i++])
          .then(address => console.log(`wrote ${i}, remaining ${writes.length - i}.`))
          .then(next)
      }
    }

    return write(writes[i++])
      .then(next)
  }).catch(onClose)
